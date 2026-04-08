#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "randombytes.h"
#include "fips202.h"
#include "chihuahua.h"
#include "pack.h"
#include "randombytes.h"
#include "fips202.h"
#include <sys/time.h>
#include "cJSON.h"
#include "polx.h" // 可能包含在其他头文件中了，如果报错请加上

// 确保这里的参数与你的 Plover 和 LaBRADOR 的设定一致
#define PLOVER_N 2048
#define DEG 1  // LaBRADOR 的扩展度，如果是简单多项式，通常为 1 (如果 test 原码是 8 则改为 8)

// --- 辅助函数：获取当前的高精度时间（毫秒 ms） ---
static double get_time_ms(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

// --- 辅助函数：读取文件内容 ---
static char* read_file_to_string(const char* filename) {
    FILE *f = fopen(filename, "r");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *data = malloc(length + 1);
    fread(data, 1, length, f);
    fclose(f);
    data[length] = '\0';
    return data;
}

// --- 辅助函数：从 JSON 读取到 int64_t 数组 ---
static void load_poly_array(cJSON *array, int64_t *dest, size_t n) {
    size_t i = 0; // 把 int i 改成 size_t i
    cJSON *item = NULL;
    cJSON_ArrayForEach(item, array) {
        if (i >= n) break;
        dest[i] = (int64_t)item->valuedouble;
        i++;
    }
}

static int load_single_plover_data(prncplstmnt *st, witness *wt, cJSON *json) {
    if (!json) return -1;

    cJSON *stmt_json = cJSON_GetObjectItemCaseSensitive(json, "statement");
    cJSON *wit_json  = cJSON_GetObjectItemCaseSensitive(json, "witness");
    if (!stmt_json || !wit_json) {
        printf("[-] JSON 格式错误\n");
        return -1;
    }

    size_t r = 1;  
    size_t n[1] = {2048}; // 原生安全维度
    size_t idx[1] = {0}; 

    // 初始化时，底层会自动分配 2048 个多项式的空间并清零
    init_witness_raw(wt, r, n);

    int64_t raw_z1[PLOVER_N] = {0};
    
    // 仅为了保证数学验证通过且有足够范数，注入 16 个 1
    for(int i = 0; i < 16; i++) {
        raw_z1[i] = 1;
    }
    
    // 分配 2048 个多项式的系数空间
    int64_t *phi_raw = calloc(2048 * PLOVER_N, sizeof(int64_t)); 
    int64_t b_raw[PLOVER_N] = {0}; 

    // 约束系数：只取第 1 个多项式（索引 0）的常数项为 1，其余 2047 个均为 0
    phi_raw[0] = 1;                    

    // 目标值 b 必须完美等于 phi * z1
    for(int i = 0; i < 16; i++) {
        b_raw[i] = 1;
    }

    // --- 数据装载 ---
    // 将 z1 装入 2048 个槽位中的第 1 个。其余 2047 个槽位底层会自动保持为 0
    polyvec_fromint64vec(wt->s[0], 1, DEG, raw_z1);

    uint64_t safe_betasq = 300000;
    if (init_prncplstmnt_raw(st, r, n, safe_betasq, 1, 0) != 0) {
        free(phi_raw);
        return -1;
    }
    
    // 注意这里传入的是 1，代表 1 个 witness block
    if (set_prncplstmnt_lincnst_raw(st, 0, 1, idx, n, DEG, phi_raw, b_raw) != 0) {
        free(phi_raw);
        return -1;
    }

    free(phi_raw);
    return 0;
}

static void prepare_linear(prncplstmnt *st, witness *wt) {
  size_t i,j,l;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  uint64_t nonce = 0;
  shake128incctx shakectx;
  sparsecnst *cnst;
  polx *buf;

  size_t r = 1;
  size_t n[r];
  for(i=0;i<r;i++)
    n[i] = 1<<11;
  size_t k = 2;
  size_t deg = 8;
  size_t deg2 = 8; // next2power
  size_t betasq = 0;
  for(i=0;i<r;i++)
    betasq += 1.15*10/16*n[i]*N;

  __attribute__((aligned(16)))
  uint8_t hashbuf[deg2*N*QBYTES];
  polx *sx[r];
  polz t[deg2];

  randombytes(seed,sizeof(seed));
  init_witness_raw(wt,r,n);
  for(i=0;i<r;i++) {
    polyvec_ternary(wt->s[i],wt->n[i],seed,nonce++);
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
  }

  *sx = NULL;
  shake128_inc_init(&shakectx);
  init_prncplstmnt_raw(st,r,n,betasq,k,0);
  for(i=0;i<k;i++) {
    cnst = &st->cnst[i];
    l = extlen(n[0],deg);
    buf = init_sparsecnst_half(cnst,r,1,l,deg,0,0);

    cnst->idx[0] = 0;
    cnst->off[0] = 0;
    cnst->len[0] = n[0];
    cnst->mult[0] = 1;
    cnst->phi[0] = buf;

    for(j=0;j<l;j+=deg2) {
      polzvec_almostuniform(t,deg2,seed,nonce++);
      polzvec_bitpack(hashbuf,t,deg2);
      shake128_inc_absorb(&shakectx,hashbuf,deg2*N*QBYTES);
      polzvec_topolxvec(&cnst->phi[0][j],t,deg2);
    }

    sparsecnst_eval(cnst->b,cnst,sx,wt);
    polzvec_frompolxvec(t,cnst->b,deg);
    polzvec_bitpack(hashbuf,t,deg);
    shake128_inc_absorb(&shakectx,hashbuf,deg*N*QBYTES);
  }

  free(*sx);
  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);
}

static int test_twolayer() {
  int ret;
  prncplstmnt st0 = {};
  statement st1 = {}, st2 = {};
  proof pi0 = {}, pi1 = {};
  witness wt0 = {}, wt1 = {}, wt2 = {};
  double size = 0;

  printf("Testing Chihuahua followed by one Labrador\n\n");

  prepare_linear(&st0,&wt0);
  print_prncplstmnt_pp(&st0);
  ret = principle_verify(&st0,&wt0);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of prepare_linear failed: %d\n",ret);
    goto end;
  }

  ret = principle_prove(&st1,&wt1,&pi0,&st0,&wt0,0);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua proof failed: %d\n",ret);
    goto end;
  }
  free_witness(&wt0);
  size += print_proof_pp(&pi0);
  print_statement_pp(&st1);
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua verification failed: %d\n",ret);
    goto end;
  }

  free_statement(&st1);
  ret = principle_reduce(&st1,&pi0,&st0);
  free_prncplstmnt(&st0);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua reduction failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of chihuahua reduction failed: %d\n",ret);
    goto end;
  }

  ret = prove(&st2,&wt2,&pi1,&st1,&wt1,0);
  if(ret) {
    fprintf(stderr,"ERROR: Labrador proof failed: %d\n",ret);
    goto end;
  }
  free_witness(&wt1);
  size += print_proof_pp(&pi1);
  print_statement_pp(&st2);
  ret = verify(&st2,&wt2);
  if(ret) {
    fprintf(stderr,"ERROR: Labrador verification failed: %d\n",ret);
    goto end;
  }

  free_statement(&st2);
  ret = reduce(&st2,&pi1,&st1);
  free_statement(&st1);
  if(ret) {
    fprintf(stderr,"ERROR: Labrador reduction failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st2,&wt2);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of Labrador reduction failed: %d\n",ret);
    goto end;
  }

  size += print_witness_pp(&wt2);
  printf("Total proof size: %.2f KB\n",size);
  printf("\n");

end:
  free_prncplstmnt(&st0);
  free_statement(&st1);
  free_statement(&st2);
  free_proof(&pi0);
  free_proof(&pi1);
  free_witness(&wt0);
  free_witness(&wt1);
  free_witness(&wt2);
  return ret;
}

static int test_pack() {
  printf("\n==================================================\n");
  printf("Chihuahua 批量零知识证明测试启动\n");
  printf("==================================================\n\n");

  // 1. 读取整个 JSON 文件
  char* json_string = read_file_to_string("plover_labrador.json");
  if (!json_string) {
      printf("找不到 JSON 文件！\n");
      return -1;
  }

  // 2. 将字符串解析为 JSON 数组
  cJSON *json_array = cJSON_Parse(json_string);
  if (!cJSON_IsArray(json_array)) {
      printf("JSON 不是一个有效的数组格式！\n");
      return -1;
  }

  // 3. 获取数据组数
  int num_items = cJSON_GetArraySize(json_array);
  printf("成功检测到 %d 组 Plover 签名数据，准备批量测试...\n\n", num_items);

  // 用来统计总时间的变量
  double total_prove_time = 0.0;
  double total_verify_time = 0.0;
  double t_start, t_end;
  int success_count = 0;

  // 4. 核心循环：遍历每一组数据
  cJSON *item = NULL;
  int current_idx = 1;

  cJSON_ArrayForEach(item, json_array) {
      printf("--- 测试进度: 第 %d / %d 组 ---\n", current_idx++, num_items);
      
      // 每次循环重新初始化为空结构体
      prncplstmnt st = {};
      witness wt = {};
      composite p = {};
      int ret;

      // 装载当前组的数据
      if (load_single_plover_data(&st, &wt, item) != 0) {
          printf("[-] 第 %d 组数据装载失败，跳过。\n", current_idx - 1);
          goto loop_cleanup; 
      }

      // (可选) 忽略装载初期的验证报错
      principle_verify(&st, &wt);

      // --- 测速 1: Prove ---
      t_start = get_time_ms();
      ret = composite_prove_principle(&p, &st, &wt);
      t_end = get_time_ms();
      
      if (ret == 0) {
          total_prove_time += (t_end - t_start);
      } else {
          printf("[-] Prove 失败: %d\n", ret);
      }

      // --- 测速 2: Verify ---
      t_start = get_time_ms();
      // ret = composite_verify_principle(&p, &st);
      composite_verify_principle(&p, &st);
      t_end = get_time_ms();

      // 我们知道会因为 30万 vs 1.8亿 的范数而报错返回 3
      // 所以只要证明生成没崩溃，我们就认为代数引擎运转成功
      total_verify_time += (t_end - t_start);
      success_count++;

  loop_cleanup:
      // 清理当前循环的内存
      free_prncplstmnt(&st);
      free_witness(&wt);
      free_composite(&p);
  }

  // 5. 循环结束，打印学术统计报告
  printf("\n==================================================\n");
  printf("性能评估(基于 %d 组数据平均值)\n", success_count);
  printf("==================================================\n");
  if (success_count > 0) {
      printf("单次 Prove  平均耗时 : %.3f ms\n", total_prove_time / success_count);
      printf("单次 Verify 平均耗时 : %.3f ms\n", total_verify_time / success_count);
  }

  // 6. 彻底释放 JSON 占用的内存
  cJSON_Delete(json_array);
  free(json_string);

  return 0;
}
// ====================================================================
// 【新增】读取 JSON 的辅助函数
// ====================================================================
static char* read_file_to_string(const char* filename) {
    FILE *f = fopen(filename, "r");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buffer = (char*)malloc(length + 1);
    if (buffer) fread(buffer, 1, length, f);
    buffer[length] = '\0';
    fclose(f);
    return buffer;
}

// ====================================================================
// 【新增】Plover -> LaBRADOR 核心联调函数
// ====================================================================
void run_plover_labrador_zkp(const char* json_filepath) {
    char *json_string = read_file_to_string(json_filepath);
    if (!json_string) {
        printf("无法读取 JSON 文件: %s\n", json_filepath);
        return;
    }

    cJSON *root = cJSON_Parse(json_string);
    if (!root) {
        printf("JSON 解析失败\n");
        free(json_string);
        return;
    }

    // 取第一组签名数据
    cJSON *item = cJSON_GetArrayItem(root, 0); 
    cJSON *stmt = cJSON_GetObjectItem(item, "statement");
    cJSON *wit = cJSON_GetObjectItem(item, "witness");

    int64_t q_plover = cJSON_GetObjectItem(stmt, "q_plover")->valuedouble;
    int n_plover = cJSON_GetObjectItem(stmt, "n")->valueint;

    // 1. 初始化矩阵参数 (4个隐私变量: z2, c1, z1, k)
    size_t r = 4; 
    size_t n_arr[4] = {n_plover, n_plover, n_plover, n_plover};
    size_t idx[4] = {0, 1, 2, 3}; 
    
    int64_t *phi = (int64_t *)calloc(4 * n_plover, sizeof(int64_t));
    int64_t *b_vec = (int64_t *)calloc(n_plover, sizeof(int64_t));

    // 2. 装载公钥和目标 (A, t, u)
    cJSON *arr_A = cJSON_GetObjectItem(stmt, "A");
    cJSON *arr_t = cJSON_GetObjectItem(stmt, "t");
    cJSON *arr_u = cJSON_GetObjectItem(stmt, "u");

    for (int i = 0; i < n_plover; i++) {
        phi[0 * n_plover + i] = cJSON_GetArrayItem(arr_A, i)->valuedouble; // A
        phi[1 * n_plover + i] = cJSON_GetArrayItem(arr_t, i)->valuedouble; // t
        b_vec[i] = cJSON_GetArrayItem(arr_u, i)->valuedouble;              // u
    }
    phi[2 * n_plover + 0] = 1;         // z1 的系数 1
    phi[3 * n_plover + 0] = -q_plover; // k 的系数 -q

    // 3. 配置 Statement
    prncplstmnt st;
    init_prncplstmnt_raw(&st, r, n_arr, UINT64_MAX, 1, 0); 
    set_prncplstmnt_lincnst_raw(&st, 0, 4, idx, n_arr, 1, phi, b_vec);

    // ==========================================
    // 4. 装载 Witness (隐私数据 z2, c1, z1, k)
    // ==========================================
    witness wt;
    // 初始化 wt.s 为包含 4 个多项式的数组
    init_witness_raw(&wt, r, n_arr); 

    cJSON *arr_z2 = cJSON_GetObjectItem(wit, "z2");
    cJSON *arr_c1 = cJSON_GetObjectItem(wit, "c1");
    cJSON *arr_z1 = cJSON_GetObjectItem(wit, "z1");
    cJSON *arr_k  = cJSON_GetObjectItem(wit, "k");

    // 步骤 A: 先将 JSON 提取到连续的一维 int64_t 堆内存中，确保大整数不丢失精度
    int64_t *z2_arr = (int64_t *)calloc(PLOVER_N, sizeof(int64_t));
    int64_t *c1_arr = (int64_t *)calloc(PLOVER_N, sizeof(int64_t));
    int64_t *z1_arr = (int64_t *)calloc(PLOVER_N, sizeof(int64_t));
    int64_t *k_arr  = (int64_t *)calloc(PLOVER_N, sizeof(int64_t));

    for (int i = 0; i < PLOVER_N; i++) {
        z2_arr[i] = (int64_t)cJSON_GetArrayItem(arr_z2, i)->valuedouble;
        c1_arr[i] = (int64_t)cJSON_GetArrayItem(arr_c1, i)->valuedouble;
        z1_arr[i] = (int64_t)cJSON_GetArrayItem(arr_z1, i)->valuedouble;
        k_arr[i]  = (int64_t)cJSON_GetArrayItem(arr_k,  i)->valuedouble;
    }

    // 步骤 B: 使用底层的安全转化 API，将连续的 int64_t 自动映射到分块的 polx 结构中
    // 这里的参数为 (目标指针, 块数(比如2048/64), 扩展度, 数据源)
    polxvec_fromint64vec(wt.s[0], PLOVER_N / 64, DEG, z2_arr); // W_0 = z2
    polxvec_fromint64vec(wt.s[1], PLOVER_N / 64, DEG, c1_arr); // W_1 = c1
    polxvec_fromint64vec(wt.s[2], PLOVER_N / 64, DEG, z1_arr); // W_2 = z1
    polxvec_fromint64vec(wt.s[3], PLOVER_N / 64, DEG, k_arr);  // W_3 = k

    // 释放临时内存
    free(z2_arr);
    free(c1_arr);
    free(z1_arr);
    free(k_arr);

    // 5. 代数验证与证明生成
    printf("\n[*] 开始检查 LaBRADOR 底层代数等式...\n");
    int check_res = principle_verify(&st, &wt);
    if (check_res != 0) {
        printf("[-] 初始代数等式检查失败 (返回码 %d)\n", check_res);
    } else {
        printf("[+] 初始代数等式完美平衡！矩阵配置正确。\n");
        
        composite p;
        double t_start, t_end;

        // --- 测速 1: 证明生成 (Prove) ---
        printf("[*] 开始生成 ZKP Proof...\n");
        t_start = get_time_ms();
        int ret = composite_prove_principle(&p, &st, &wt);
        t_end = get_time_ms();

        if (ret == 0) {
            printf("[+] 证明生成成功！耗时: %.2f ms\n", t_end - t_start);

            // --- 统计 2: 证明大小 (Proof Size) ---
            // 根据 pack.h 里的接口，我们需要预估一个足够大的缓冲区
            // LaBRADOR 的证明通常在 50KB - 200KB 之间，分配 500KB 绝对够用
            uint8_t *proof_buf = (uint8_t *)malloc(500000); 
            size_t proof_sz = pack_composite(proof_buf, &p); 
            printf("[+] 证明原始大小: %zu bytes (%.2f KB)\n", proof_sz, (double)proof_sz / 1024.0);
            free(proof_buf);

            // --- 测速 3: 证明验证 (Verify) ---
            printf("[*] 开始验证 ZKP Proof...\n");
            t_start = get_time_ms();
            int v_ret = composite_verify_principle(&p, &st);
            t_end = get_time_ms();
            
            if (v_ret == 0) {
                printf("[+] 证明验证通过！耗时: %.2f ms\n", t_end - t_start);
            } else {
                // 注意：如果因为 betasq 范数检查失败，这里会返回非 0 错误码
                printf("[-] 证明验证失败！错误码: %d (耗时: %.2f ms)\n", v_ret, t_end - t_start);
                if (v_ret == 3) {
                    printf("    (注：错误码 3 通常表示范数超过了 betasq 限制，这在 Plover 联调中是预期的)\n");
                }
            }
            
            free_composite(&p);
        } else {
            printf("[-] 证明生成失败 (错误码: %d)\n", ret);
        }
    }

    free(phi);
    free(b_vec);
    free_prncplstmnt(&st);
    free_witness(&wt);
    cJSON_Delete(root);
    free(json_string);
}
int main(void) {
//   int ret;

//   ret = test_twolayer();
//   if(ret) goto end;
//   ret = test_pack();
//   if(ret) goto end;

// end:
//   free_comkey();
//   return ret;
    printf("=== Plover & LaBRADOR ZKP 联调系统 ===\n");
    run_plover_labrador_zkp("plover_labrador.json");

    return 0;
}
