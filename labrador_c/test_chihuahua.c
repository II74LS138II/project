#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#include "randombytes.h"
#include "fips202.h"
#include "chihuahua.h"
#include "pack.h"
#include "cJSON.h"
#include "polx.h"

// 确保与 Plover 参数严格对齐
#define PLOVER_N 2048
#define DEG 1  

// 函数原型声明
void run_plover_labrador_zkp(const char* json_filepath);
static double get_time_ms(void);
static char* read_file_to_string(const char* filename);

// --- 1. 高精度计时器 ---
static double get_time_ms(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

// --- 2. 安全读取文件 ---
static char* read_file_to_string(const char* filename) {
    FILE *f = fopen(filename, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *data = malloc(length + 1);
    if (data) {
        size_t r = fread(data, 1, length, f);
        data[r] = '\0';
    }
    fclose(f);
    return data;
}

// --- 3. 核心联调逻辑 ---
void run_plover_labrador_zkp(const char* json_filepath) {
    char *json_string = read_file_to_string(json_filepath);
    if (!json_string) {
        printf("[-] 无法读取 JSON 文件: %s\n", json_filepath);
        return;
    }

    cJSON *root = cJSON_Parse(json_string);
    if (!root) {
        printf("[-] JSON 解析失败\n");
        free(json_string);
        return;
    }

    // 取第一组签名数据进行验证
    cJSON *item = cJSON_GetArrayItem(root, 0); 
    cJSON *stmt_json = cJSON_GetObjectItem(item, "statement");
    cJSON *wit_json = cJSON_GetObjectItem(item, "witness");

    int64_t q_plover = (int64_t)cJSON_GetObjectItem(stmt_json, "q_plover")->valuedouble;
    
    // 1. 初始化矩阵参数 (4个隐私变量: z2, c1, z1, k)
    size_t r = 4; 
    size_t n_arr[4] = {PLOVER_N, PLOVER_N, PLOVER_N, PLOVER_N};
    size_t idx[4] = {0, 1, 2, 3}; 
    
    int64_t *phi = (int64_t *)calloc(4 * PLOVER_N, sizeof(int64_t));
    int64_t *b_vec = (int64_t *)calloc(PLOVER_N, sizeof(int64_t));

    // 2. 装载公钥 A, t 和目标 u
    cJSON *arr_A = cJSON_GetObjectItem(stmt_json, "A");
    cJSON *arr_t = cJSON_GetObjectItem(stmt_json, "t");
    cJSON *arr_u = cJSON_GetObjectItem(stmt_json, "u");

    for (int i = 0; i < PLOVER_N; i++) {
        phi[0 * PLOVER_N + i] = (int64_t)cJSON_GetArrayItem(arr_A, i)->valuedouble; // A
        phi[1 * PLOVER_N + i] = (int64_t)cJSON_GetArrayItem(arr_t, i)->valuedouble; // t
        b_vec[i] = (int64_t)cJSON_GetArrayItem(arr_u, i)->valuedouble;              // u
    }
    phi[2 * PLOVER_N + 0] = 1;         // z1 的系数多项式为 1
    phi[3 * PLOVER_N + 0] = -q_plover; // k 的系数多项式为 -q

    // 3. 配置 Statement (放开范数限制)
    prncplstmnt st;
    uint64_t infinite_beta = (uint64_t)-1; // 解决 Total witness norm too big
    if (init_prncplstmnt_raw(&st, r, n_arr, infinite_beta, 1, 0) != 0) {
        printf("[-] Statement 初始化失败\n");
        goto cleanup_json;
    }
    set_prncplstmnt_lincnst_raw(&st, 0, 4, idx, n_arr, DEG, phi, b_vec);

    // 4. 装载 Witness (z2, c1, z1, k)
    witness wt;
    init_witness_raw(&wt, r, n_arr); 

    int64_t *tmp = malloc(PLOVER_N * sizeof(int64_t));
    const char *wit_names[4] = {"z2", "c1", "z1", "k"};
    for(int j=0; j<4; j++) {
        cJSON *arr = cJSON_GetObjectItem(wit_json, wit_names[j]);
        for(int i=0; i<PLOVER_N; i++) tmp[i] = (int64_t)cJSON_GetArrayItem(arr, i)->valuedouble;
        polxvec_fromint64vec((polx *)wt.s[j], PLOVER_N/64, DEG, tmp);
    }
    free(tmp);

    // 5. 执行证明与验证
    printf("\n[*] 开始检查 LaBRADOR 底层代数等式...\n");
    if (principle_verify(&st, &wt) != 0) {
        printf("[-] 初始代数等式检查失败！请检查 Python 端的 k 提取逻辑。\n");
    } else {
        printf("[+] 初始代数等式完美平衡！\n");
        
        composite p;
        double t_start, t_end;

        printf("[*] 开始生成 ZKP Proof...\n");
        t_start = get_time_ms();
        int ret = composite_prove_principle(&p, &st, &wt);
        t_end = get_time_ms();

        if (ret == 0) {
            printf("[+] 证明生成成功！耗时: %.2f ms\n", t_end - t_start);
            printf("[+] 证明压缩大小: %.2f KB\n", p.size);

            printf("[*] 开始验证 ZKP Proof...\n");
            t_start = get_time_ms();
            int v_ret = composite_verify_principle(&p, &st);
            t_end = get_time_ms();
            printf("[%s] 证明验证结果: %d (耗时: %.2f ms)\n", v_ret==0?"+":"-", v_ret, t_end - t_start);
            
            free_composite(&p);
        } else {
            printf("[-] 证明生成失败: %d\n", ret);
        }
    }

    // 资源释放
    free_prncplstmnt(&st);
    free_witness(&wt);
    free(phi);
    free(b_vec);

cleanup_json:
    cJSON_Delete(root);
    free(json_string);
}

// --- 4. Main 入口 ---
int main(void) {
    // 初始化 LaBRADOR 随机源
    uint8_t entropy[48];
    
    // 1. 调用你 randombytes.h 中声明的唯一函数
    randombytes(entropy, 48);

    printf("=== Plover & LaBRADOR 跨语言 ZKP 联调系统 ===\n");
    
    // 处理 Python 导出的 JSON
    run_plover_labrador_zkp("plover_export.json");

    return 0;
}