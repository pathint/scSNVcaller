#include <stdio.h>
#include <string.h>
#include <htslib/sam.h>
#include <omp.h>
#define MAX_CELLS 10000
	
int main(int argc, char *argv[]){
	
	// 0、如果terminal输入的参数不够，提示用户使用方法
	if (argc  < 4){
		printf("Usage: %s <input.bam> <Sample-Barcode.File> <output.bam>",argv[0]);
	}
	
    char targetCells[MAX_CELLS][20];            // 二维字符数组，用于存放读入文件的细胞ID，第一维表示最大细胞数目，第二维表示每个细胞ID的长度 
    char targetCell[20];                        // 用于存取当前行要获取的细胞ID
	
	// 1、读入样本-细胞文件并将其转为字符串数组	
    char *filename = argv[2];                   // 将程序运行时传入的第3个参数存储到变量filename中
    FILE *sample_barcode;                       // 定义了一个指向FILE类型的指针变量，
    sample_barcode = fopen(filename,"r");       // 用fopen函数 打开名为‘filename’的文件，并将其赋值给 sample_barcode，“r”表示以只读方式打开文件
	
    if (sample_barcode == NULL){    // 样本细胞-文件打开失败，退出运行
        printf("Failed to open file %s\n",filename);
        return 1;
    }   
    char line[200];                             // 定义一个名为‘line’的字符数组，用于存储每一行读入的文件内容，最大容纳字符数为200
    int i = 0;                                  // 定义整型变量i，用于存储当前已经读取的行数
    while (fgets(line, 200,sample_barcode) != NULL  ){  // 用于读取文件中的每一行数据，fgets用于从文件中读取一行数据，第一个参数是一个字符数组，用于存储读取的数据，第二个参数是最大读取的字符数，第三个参数是要读取的文件指针。 如果读取成功，则会返回读取的字符串指针，否则返回NULL。while循环会不断读取文件中的每一行数据，知道读取到文件末尾
        sscanf(line,"%s",&targetCell);          // 将每行读到的第一列赋值给targetCell   
        strcpy(targetCells[i],targetCell);      // 将targetCell变量中的字符串复制到字符串数组targetCells中的第i个位置
        i++;                                    // i自增1，表示已经读取了一行数据
    }
	
	// 2、bam文件的操作
	samFile *fp_in = sam_open(argv[1],"r");		// 打开bam/sam文件
	bam_hdr_t *bhdr = sam_hdr_read(fp_in);			// 读取bam/sam文件头信息
	bam1_t *b = bam_init1();					// 初始化一个bam1_t1结构体，用于存储读入的一个bam记录
	kstring_t str = KS_INITIALIZE;              // 用于存储整个bam文件记录	

	char *output_fname = argv[3];
    samFile *fp_out = sam_open(output_fname,"wb");
    if (sam_hdr_write(fp_out, bhdr) < 0) {
        fprintf(stderr, "Error: failed to write header to file %s\n", output_fname);
        exit(1);
    }   

	int ret;
	int k = 0;									// 当前查找的是样本-细胞文件里的第k个细胞 
	
	strcpy(targetCell,targetCells[k]);			// 刚开始时，待提取的第一个细胞
	

	int x = 1;
	while (sam_read1(fp_in,bhdr,b) >= 0){ 
    	uint8_t *cell_index = bam_aux_get(b,"CB");
		const char *cell;   
        if (cell_index){                        // 如果是真实有效的 Cell Barcode，则进行下一步分析
            cell = bam_aux2Z(cell_index);       // 获得细胞的Cell Barcode
	
			if (strcmp(cell,targetCell) < 0 && k < i){
				continue;
			}else if (strcmp(cell,targetCell) == 0 && k < i ){
				if (sam_write1(fp_out,bhdr,b) < 0 ){  
            		fprintf(stderr,"Error: failed to write alignment to output file %s\n",output_fname);
            		exit(1);
				}
			}else if (strcmp(cell,targetCell) >0 && k < i ){
				k++;
				strcpy(targetCell,targetCells[k]);
			}
		}
		if (k >= i){
			break;
		}
	}
    bam_destroy1(b);
    bam_hdr_destroy(bhdr);
    sam_close(fp_in);
    sam_close(fp_out);
	fclose(sample_barcode);
	return 0;
}
