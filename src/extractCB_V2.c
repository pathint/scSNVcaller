#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <string.h>
#define MAX_CELLS 1000 
#define MAX_LINE_LENGTH 100	
	

int main(int argc, char *argv[]){
	
	// 0、如果terminal输入的参数不够，提示用户使用方法
	if (argc < 4){
		printf("Usage: %s <bamfile> <vcffile> <Sample>\n", argv[0]);
		return 1;
	}
	
	// 1、vcf文件的操作
	htsFile *vcf = bcf_open(argv[2],"r");   
	bcf_hdr_t *vhdr = bcf_hdr_read(vcf); 
	bcf1_t *rec = bcf_init1();			 
	char location[200];							// 声明一个变量，用于存放vcf中获得的”染色体-位点“，作为bam文件处理时的输入，并为其分配内存空间
	
	// 2、bam文件的操作
	samFile *fp = sam_open(argv[1],"r");
	bam_hdr_t *bhdr = sam_hdr_read(fp);	
	bam1_t *b = bam_init1();		
	
		// bam文件遍历	
	hts_idx_t *idx;	
	hts_itr_t *iter;
	
	int start;									// 存放当前对比序列的起始位点-1	
		// 样本名：第4个输入参数
	char sample[100];
	strcpy(sample, argv[3]);
	
	
	// 3、输出文件
	char refFile[200];													// 存放没有突变的细胞的文件
	char altFile[200];													// 存放突变的细胞的文件
	char normalCells[MAX_CELLS * MAX_LINE_LENGTH];						// 用于正常存储细胞的大字符串
	char tumorCells[MAX_CELLS * MAX_LINE_LENGTH];						// 用于存储肿瘤细胞的大字符串

	int nc;																// 正常细胞计数
	int tc;																// 肿瘤细胞计数
	
	while (bcf_read(vcf,vhdr,rec) == 0){			// 遍历vcf文件	
		normalCells[0] = '\0';                                          // 将大字符串初始化为空字符串
	    tumorCells[0] = '\0';                                           // 将大字符串初始化为空字符串

		bcf_unpack(rec,BCF_UN_ALL);							
		const char *chrom = bcf_hdr_id2name(vhdr,rec->rid);	
		int pos = rec->pos + 1;								
		
		snprintf(location,sizeof(location),"%s:%d-%d",chrom,pos,pos);	// 待查询位点
		
		idx = sam_index_load(fp,argv[1]);								// 加载索引信息
		iter = sam_itr_querys(idx,bhdr,location);						// 从索引信息获得执行的比对序列
				
		char *ref = rec->d.allele[0];						
		char *alt = rec->d.allele[1];						
		kstring_t str = KS_INITIALIZE;	
		
		// 将变量名和文件扩展名组合成文件名
		snprintf(refFile,sizeof(refFile),"%s_%s_%d_ref_%s.txt",sample,chrom,pos,ref);	
		snprintf(altFile,sizeof(altFile),"%s_%s_%d_alt_%s.txt",sample,chrom,pos,alt);	

		// 打开输出文件，如果没有输出文件，就创建一个输出文件
		FILE *refOut = fopen(refFile,"w");
		FILE *altOut = fopen(altFile,"w");
		
			
		while (sam_itr_next(fp,iter,b) >= 0){		// 遍历 待查询 比对序列
	
			uint8_t *cell_index = bam_aux_get(b,"CB");
		
			char *cell;									
			if (cell_index){						// 如果是真实有效的 Cell Barcode，则进行下一步分析
				cell = bam_aux2Z(cell_index);							// 获取有效的比对序列对应的细胞的Cell Barcode
				
				// CIGAR 字符串的处理
				uint32_t *cigar = bam_get_cigar(b);						// 获取CIGAR字符串
				int n_cigar = b->core.n_cigar;							// 获取当前CIGAR字符串包含的操作个数
				uint8_t *seq = bam_get_seq(b);							// 获取序列
				int32_t seqlen =  bam_cigar2qlen(n_cigar,cigar);		// 获取比对序列的长度
				start = b->core.pos;									// 可变的起始位点-1
				int k ;													// 当前CIGAR的第k个操作
	
				int j = 0;												// 用于存储找到的基因是比对序列上的第几个碱基
				int offset;
				int base_num;											// 用于存储当前碱基在字符串“=ACMGRSVTWYHKDBN”的序号
				cell = bam_aux2Z(cell_index);							// 获得当前比对序列对应的Cell Barcode

				for(k=0; k < n_cigar; k++){			// 遍历CIGAR字符串中的每个操作
					int oplen = bam_cigar_oplen(cigar[k]);				// 用于记录CIGAR中当前第k个操作的操作长度
					int op = bam_cigar_op(cigar[k]);					// 用于记录CIGAR中当前第k个操作的操作类型（数值）
					char opchr = bam_cigar_opchr(cigar[k]);				// 用于记录CIGAR中当前第k个操作的操作类型
					
					if (start + oplen < pos){		// 如果 “可变起始位点”+”操作长度“ < 查询位点  即 不覆盖的情况
						if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CINS  || op == BAM_CDIFF || op == BAM_CSOFT_CLIP){  // M、=、I、S、X
							j += oplen; 
						}
						if (op != BAM_CINS ){
							start += oplen; 							// 可变的起始位点 变化
							if (op == BAM_CSOFT_CLIP && k == 0){
								start -= oplen;
							}
						}
					}else{							// 覆盖到 查询位点
						if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CINS  || op == BAM_CDIFF || op == BAM_CSOFT_CLIP){  // M、=、I、S、X
							if((op == BAM_CSOFT_CLIP || op == BAM_CINS ) && k == 0){
								j += oplen;
								continue;	
							}
			
							offset = pos - start;
							j += offset - 1;
							base_num = bam_seqi(seq,j);					// 获取当前碱基在字符串“=ACMGRSVTWYHKDBN”的序号
							char base = "=ACMGRSVTWYHKDBN"[base_num]; 
							if (base == *ref){		// 如果当前碱基与参考基因一致
								if (strstr(normalCells,cell) == NULL){
									strcat(normalCells,cell);
									fprintf(refOut,"%s\n",cell);
								}
							}else if (base == *alt){
								if (strstr(tumorCells,cell) == NULL){
									strcat(tumorCells,cell);
									fprintf(altOut,"%s\n",cell);
								}
							}
						}	
						break;											// 在覆盖到查询位点之后，就不再进行下一个CIGAR字符的循环了
					}
					
				}
			}
				
		}
		
	}
	// 关闭已打开的文件
		// （1） vcf文件关闭
	bcf_destroy1(rec);
	bcf_hdr_destroy(vhdr);
	bcf_close(vcf);
		//  (2)  bam文件关闭
	bam_destroy1(b);	
	bam_hdr_destroy(bhdr);
	hts_itr_destroy(iter);
	hts_idx_destroy(idx);
	sam_close(fp);

	return 0;
			
}

