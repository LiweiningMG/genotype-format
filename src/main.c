#include <argp.h>
#include <stdio.h>
#include <stdlib.h>

#include "common_fun.h"
#include "ped_format.h"

#define MAX_SAMPLE_NUM 5000
#define MAX_ID_LEN 10

int main(int argc, char *argv[]) {
  /* ----------- 命令行参数处理 ----------- */
  struct arguments para;
  //   int header = 1;
  //   int rowname = 1;

  /* 可选参数默认值 */
  para.out = (char *)"plink.ped";

  /* 解析命令行参数  */
  para.require_para = 0; /* 必要参数个数 */
  argp_parse(&argp, argc, argv, ARGP_IN_ORDER, 0, &para);

  /* 基因型文件中个体数 */
  para.nMrk = countLines(para.geno) - 1;

  char *file_char = read_whole(para.geno);

  // id提取
  char *all_ptr = NULL;
  char *row_ptr = NULL;
  char *value = NULL;
  char *row = NULL;
  char **id = (char **)calloc(MAX_SAMPLE_NUM * MAX_ID_LEN, sizeof(char *));
  int *id_len = (int *)calloc(MAX_SAMPLE_NUM, sizeof(int));
  row = strtok_r(file_char, "\n", &all_ptr);
  value = strtok_r(row, "\t", &row_ptr);

  // 有行名
  if (1) {
    value = strtok_r(NULL, "\t", &row_ptr);
  }

  // 保存个体ID
  para.nInd = 0;
  while (value != NULL) {
    if (para.nInd >= MAX_SAMPLE_NUM) {
      fprintf(stderr, "Too many samples, max %d\n", MAX_SAMPLE_NUM);
      return 1;
    }

    id[para.nInd] = (char *)calloc(MAX_ID_LEN, sizeof(char));
    strcpy(id[para.nInd], value);
    id_len[para.nInd] = strlen(id[para.nInd]);

    para.nInd++;
    value = strtok_r(NULL, "\t", &row_ptr);
  }

  // 准备ped文件前六列
  int id_len_all = sumi(id_len, para.nInd, 1);
  int *gt_pos = (int *)calloc(para.nInd, sizeof(int));
  int ped_len = id_len_all * 2 + para.nInd * (para.nMrk * 4 + 11);
  char *ped = (char *)calloc(ped_len, sizeof(char));
  for (size_t i = 0; i < ped_len; i++) {
    // 初始化指针数组
    ped[i] = 32; // ASCII code for spaces
  }
  int pos = 0;
  for (size_t i = 0; i < para.nInd; i++) {
    strcpy(ped + pos, id[i]);
    ped[pos + id_len[i]] = 32;
    strcpy(ped + pos + id_len[i] + 1, id[i]);
    ped[pos + id_len[i] * 2 + 1] = 32;

    free(id[i]);

    pos += id_len[i] * 2 + 1;
    ped[pos + 1] = 48; // ASCII code for 0
    ped[pos + 3] = 48;
    ped[pos + 5] = 48;
    ped[pos + 7] = 45; // ASCII code for -
    ped[pos + 8] = 57; // ASCII code for 9
    pos += 8 + 1;

    // 每个个体snp基因型起始的位置(' A B'的第一个空格位置)
    gt_pos[i] = pos;

    pos += para.nMrk * 4;
    ped[pos] = 10; // ASCII code for \n
    pos += 1;
  }
  free(id);

  // 准备ped中基因型部分
  for (size_t i = 0; i < para.nMrk; i++) {
    row = strtok_r(NULL, "\n", &all_ptr);
    // snp id
    value = strtok_r(row, "\t", &row_ptr);
    for (size_t j = 0; j < para.nInd; j++) {
      value = strtok_r(NULL, "\t", &row_ptr); // 基因型

      if (strcmp(value, "NA") == 0) {
        /* value[0] = 48; */
        value[0] = 48;
        value[1] = 48;
      }

      ped[gt_pos[j] + i * 4 + 1] = value[0];
      ped[gt_pos[j] + i * 4 + 3] = value[1];
    }
  }

  FILE *fp_out = fopen(para.out, "w");
  if (fp_out == NULL) {
    fprintf(stderr, "Failed to open file %s\n", para.out);
    return 1;
  }

  fprintf(fp_out, ped);
  fclose(fp_out);
  free(id_len);
  free(gt_pos);
  free(ped);

  return 0;
}
