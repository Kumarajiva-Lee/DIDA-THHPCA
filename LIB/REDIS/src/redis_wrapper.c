#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hircluster.h>

#include "logger.h"

// 此处保存错误信息
char err_str[256]="\0";

struct obsF{
    float lon, lat, lev;
    int id;
};

struct obsD{
    double lon, lat, lev;
    int id;
};

// 设置 socket 超时
int set_redis_socket_timeout(redisClusterContext *cc, int sec) {
    struct timeval timeout = {sec, 0};
    return redisClusterSetOptionTimeout(cc, timeout);
}

// 以下为错误信息处理
void err_clear() {
    sprintf(err_str, "");
    return;
}

void err_print(char *fmt, ...) {
    va_list arg_ptr;
    va_start(arg_ptr, fmt);
    vsprintf(err_str, fmt, arg_ptr); 
    va_end(arg_ptr);
}

void err_append(char *fmt, ...) {
    int l = strlen(err_str);
    va_list arg_ptr;
    va_start(arg_ptr, fmt);
    vsprintf(err_str + l, fmt, arg_ptr); 
    va_end(arg_ptr);
}

void get_err_string(char *err) {
    sprintf(err, "%s", err_str);
}

void err_to_stdout() {
    printf("Error: %s", err_str);
}

// 针对不同返回类型进行不同处理
int reply_returnf(redisClusterContext *cc, redisReply *reply, float *buf, int *tag);
int reply_returnd(redisClusterContext *cc, redisReply *reply, double *buf, int *tag);
int reply_returni(redisClusterContext *cc, redisReply *reply, int *buf, int *tag);
int reply_returns(redisClusterContext *cc, redisReply *reply, char *buf, int *tag);
int reply_check(redisClusterContext *cc, redisReply *reply);

// 获取 nodes 中第 index 个集群节点
void get_ith_node(int index, char *nodes, char *node) {
    int len = strlen(nodes);
    int i;

    int nodes_cnt = 0;
    for (i = 0; i < len; i++)
        if (nodes[i] == ',')
            nodes_cnt++;

    index = index % nodes_cnt;

    i = 0;
    int split_cnt = 0;
    int s, e;
    while (split_cnt < index) {
        if (nodes[i] == ',')
            split_cnt++;
        i++;
    }
    s = i;
    while (split_cnt < index + 1 && i <= len) {
        if (nodes[i] == ',')
            split_cnt++;
        i++;
    }
    e = i - 2;
    strncpy(node, nodes + s, e - s + 1);
    node[e - s + 1] = '\0';
}

// 集群均衡连接，不同进程分散连接
redisClusterContext* redisClusterConnectBalanced(char *nodes, int rank) {
    char node[50];
    get_ith_node(rank, nodes, node);
    return redisClusterConnect(node, 0);
}

// 将 cc 缓冲区中所有命令发送
int redis_sendall(redisClusterContext *cc) {
    int ret = redisClusterSendAll(cc); 
    if (ret == REDIS_ERR) {
        err_print("EXECUTE_COMMAND_ERROR:send all error:%s", cc->errstr);
		return -1;
    }
	return 0;
}

// 测试重传机制时使用
int redis_sendall_robust_test(redisClusterContext *cc) {
    int ret = redisClusterSendAllRobust(cc); 
    if (ret == REDIS_ERR) {
        err_print("EXECUTE_COMMAND_ERROR:send all error:%s", cc->errstr);
		return -1;
    }
	return 0;
}

// 裁剪字符串
void trim_string(char *str)
{
    char *start, *end;
    int len = strlen(str);
 
    if(str[len-1] == '\n') {
        len--;
        str[len] = 0;
    }

    start = str;
    end = str + len -1;
    while(*start && isspace(*start))
        start++;
    while(*end && isspace(*end))
        *end-- = 0;
    strcpy(str, start);
}

// 与集群建立连接，会连接同一个节点
redisClusterContext *setupConnection(char *redis_address)
{
    redisReply   *reply;

    redisClusterContext *cc;
    cc = redisClusterContextInit();
	cc->flags |= REDIS_BLOCK;
    redisClusterSetOptionAddNodes(cc, redis_address);
    redisClusterSetOptionRouteUseSlots(cc);
    redisClusterConnect2(cc);
    if (cc == NULL || cc->err) {
        if (cc) {
            printf("Connection error: %s\n", cc->errstr);
            err_print("Connection error: %s\n", cc->errstr);
            redisClusterFree(cc);
        }
        else {
            printf("Connection error: can't allocate redis context\n");
            err_print("Connection error: can't allocate redis context\n");
        }
		return NULL;
    }

    return cc;
}

// 字符串类型 key value 
int redis_set(redisClusterContext *cc, char *key, char *value)
{
    redisReply *reply = redisClusterCommand(cc, "SET %s %s", key, value);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
	return errno;
}

// 二进制类型 key value 
int redis_setb(redisClusterContext *cc, char *key, void *value, int * len)
{
    redisReply *reply = redisClusterCommand(cc, "SET %s %b", key, value, *len);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
	return errno;
}

// 获取 key 对应的 value
int redis_get(redisClusterContext *cc, char *key, void *value)
{
    redisReply *reply = redisClusterCommand(cc, "GET %s", key);
	if (reply == NULL)
		return -1;
	int errno = 0;
    switch(reply->type) {
        case 1 :
        case 3 :
            memcpy(value, reply->str, reply->len);
            break;
        case 4:
            printf("return null\n");
			errno = -2;
            break;
        case 5:
            printf("return status: %s\n", reply->str);
			errno = -2;
            break;
        case 6:
            printf("return error: %s\n", reply->str);
			errno = -2;
            break;
        default:
            printf("return no match (error)!\n");
			errno = -2;
    }
    freeReplyObject(reply);
	return errno;
}

// 执行 cmd
int redis_cmd(redisClusterContext *cc, char *cmd)
{
    redisReply *reply = redisClusterCommand(cc, cmd);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
	return errno;
}

// 执行返回值为 double 的 cmd
int redis_cmdreturnd(redisClusterContext *cc, char *cmd, double *buf)
{
    int i;
    int mm=0;
    redisReply *reply = redisClusterCommand(cc, cmd);
    int errno = reply_returnd(cc, reply, buf, &mm);
    freeReplyObject(reply);
	return errno;
}

// 执行返回值为 float 的 cmd
int redis_cmdreturnf(redisClusterContext *cc, char *cmd, float *buf)
{
    int i;
    int mm=0;
    redisReply *reply = redisClusterCommand(cc, cmd);
    int errno = reply_returnf(cc, reply, buf, &mm);
    freeReplyObject(reply);
	return errno;
}


// 执行返回值为 inf 的 cmd
int redis_cmdreturni(redisClusterContext *cc, char *cmd, int *buf)
{
    int i;
    int mm=0;
    redisReply *reply = redisClusterCommand(cc, cmd);
    int errno = reply_returni(cc, reply, buf, &mm);
    freeReplyObject(reply);
	return errno;
}

void redistestset(redisClusterContext *c)
{
    printf("attempting to get wxsc\n");
    redisReply *reply = redisClusterCommand(c, "get wxsc");
    printf("get value: %s\n", reply->str);
}

// 是否使用二进制接口
#ifdef BINARY

// 以下用于设定键值的长度限制
#define MAX_KEY_LEN 100
#define MAX_FIELD_LEN 20

// hash 表设置单个int value
int RedisHseti(redisClusterContext *cc, char *hash, char *key, int *value)
{
    int argc = 4;
    char* argv[4];
    argv[0] = "hset";
    argv[1] = hash;
    argv[2] = key;
    argv[3] = (char *)value;

    size_t argvlen[4];
    argvlen[0] = 4;
    argvlen[1] = strlen(hash);
    argvlen[2] = strlen(key);
    argvlen[3] = sizeof(int);

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
	return errno;
}


// hash 表获取单个int value
int RedisHgeti(redisClusterContext *cc, char *hash, char *key, int *value)
{
    int argc = 3;
    char* argv[3];
    argv[0] = "hget";
    argv[1] = hash;
    argv[2] = key;

    size_t argvlen[3];
    argvlen[0] = 4;
    argvlen[1] = strlen(hash);
    argvlen[2] = strlen(key);

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
	int errno = reply_returni(cc, reply, value, 0);
    freeReplyObject(reply);
	return errno;
}

// hash 表设置单个float value
int RedisHsetf(redisClusterContext *cc, char *hash, char *key, float *value)
{
    int argc = 4;
    char* argv[4];
    argv[0] = "hset";
    argv[1] = hash;
    argv[2] = key;
    argv[3] = (char *)value;

    size_t argvlen[4];
    argvlen[0] = 4;
    argvlen[1] = strlen(hash);
    argvlen[2] = strlen(key);
    argvlen[3] = sizeof(float);

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
	int errno = reply_check(cc, reply);
    freeReplyObject(reply);
	return errno;
}

// hash 表获取单个float value
int RedisHgetf(redisClusterContext *cc, char *hash, char *key, float *value)
{
    int argc = 3;
    char* argv[3];
    argv[0] = "hget";
    argv[1] = hash;
    argv[2] = key;

    size_t argvlen[3];
    argvlen[0] = 4;
    argvlen[1] = strlen(hash);
    argvlen[2] = strlen(key);

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
	int errno = reply_returnf(cc, reply, value, 0);
    freeReplyObject(reply);
	return errno;
}

// hash 表设置单个double value
int RedisHsetd(redisClusterContext *cc, char *hash, char *key, double *value)
{
    int argc = 4;
    char* argv[4];
    argv[0] = "hset";
    argv[1] = hash;
    argv[2] = key;
    argv[3] = (char *)value;

    size_t argvlen[4];
    argvlen[0] = 4;
    argvlen[1] = strlen(hash);
    argvlen[2] = strlen(key);
    argvlen[3] = sizeof(double);

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
	return errno;
}

// hash 表获取单个double value
int RedisHgetd(redisClusterContext *cc, char *hash, char *key, double *value)
{
    int argc = 3;
    char* argv[3];
    argv[0] = "hget";
    argv[1] = hash;
    argv[2] = key;

    size_t argvlen[3];
    argvlen[0] = 4;
    argvlen[1] = strlen(hash);
    argvlen[2] = strlen(key);

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
	int errno = reply_returnd(cc, reply, value, 0);
    freeReplyObject(reply);
	return errno;
}

// hash 表设置单个string value
int RedisHsets(redisClusterContext *cc, char *hash, char *key, char *value)
{
    redisReply *reply = redisClusterCommand(cc, "HSET %s %s %s", hash, key, value);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
	return errno;
}

// hash 表获取单个string value
int RedisHgets(redisClusterContext *cc, char *hash, char *key, char *value, int *len1, int *len2, int *len3)
{
    redisReply *reply = redisClusterCommand(cc, "HGET %s %s", hash, key);
	int errno = reply_returns(cc, reply, value, len3);
    freeReplyObject(reply);
	return errno;
}

// hash 表设置一维int物理场
int redis_hmsetf1d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, float *buffer)
{
    int i, j, v;
    int n = *nn;
    int kv_cnt = ((*ite - *its) / *istride + 1) * n;

    int argc = kv_cnt * 2 + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    /*string split*/
    char var[n][32];
    char* token = strtok(varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for (i = *its; i <= *ite; i += *istride) {
        for (v = 0; v < n; v++) {
            sprintf(keys + offset, "%d:%s", i, var[v]);
            argv[pos * 2 + 2] = keys + offset;
            argvlen[pos * 2 + 2] = strlen(keys + offset);
            offset += argvlen[pos * 2 + 2];
            argv[pos * 2 + 3] = (char *)(buffer + ((i - *its) * n + v));
            argvlen[pos * 2 + 3] = sizeof(float);
            pos++;
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    int errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
}


// hash 表获取一维int物理场
int redis_hmgetf1d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, float *buffer)
{
    int i, j, v;
    int n = *nn;
    int kv_cnt = ((*ite - *its) / *istride + 1) * n;

    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    /*string split*/
    char var[n][32];
    char* token = strtok(varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for (i = *its; i <= *ite; i += *istride) {
        for (v = 0; v < n; v++) {
            sprintf(keys + offset, "%d:%s", i, var[v]);
            argv[pos + 2] = keys + offset;
            argvlen[pos + 2] = strlen(keys + offset);
            offset += argvlen[pos + 2];
            pos++;
        }
    }

    int mm=0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    int errno = reply_returnf(cc, reply, buffer, &mm); 

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
} 


// hash 表设置二维float物理场
int redis_hmsetf2d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, float *buffer)
{
    int i, j, v;
    int n = *nn;
	int errno = 0;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v =0; v<n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
				logger_info("var", "%d:0:%s", j,var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(float);
                pos++;
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
			errno = -1;
			goto err;
    	}
    }


    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("redisClusterGetReply Error!\n");
        	freeReplyObject(reply);
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
			errno = -1;
			goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
		if (errno) 
			goto err;
    }
err:
    free(argv);
    free(argvlen);
    free(keys);
	return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(float);
                pos++;
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
#endif
}

// hash 表获取二维float物理场
int redis_hmgetf2d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, float *buffer)
{
    int i, j, v;
	int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];
  
    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
	
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
		int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
			errno = -1;
			goto err;
    	}
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
			errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
			freeReplyObject(reply);
			goto err;
        }
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
		if (errno)
			goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
	return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v=0; v<n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnf(cc, reply, buffer, &mm);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
#endif
}

// hash 表设置三维float物理场
int redis_hmsetf(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
	int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(float);
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
			errno = -1;
			goto err;
    	}
    }

    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
			errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
			freeReplyObject(reply);
			goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
		if (errno) 
			goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
	return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(float);
                    pos++;
                }
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
#endif
}

// hash 表获取三维float物理场
int redis_hmgetf(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
	int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
		offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
			errno = -1;
			goto err;
    	}
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
			errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
			freeReplyObject(reply);
			goto err;
        }
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
		if (errno)
			goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
	return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnf(cc, reply, buffer, &mm);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
#endif
}

// hash 表设置一维double物理场
int redis_hmsetd1d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, double *buffer)
{
    int i, j, v;
	int errno = 0;
    int n = *nn;
    int kv_cnt = ((*ite - *its) / *istride + 1) * n;

    int argc = kv_cnt * 2 + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    /*string split*/
    char var[n][32];
    char* token = strtok(varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for (i = *its; i <= *ite; i += *istride) {
        for (v = 0; v < n; v++) {
            sprintf(keys + offset, "%d:%s", i, var[v]);
            argv[pos * 2 + 2] = keys + offset;
            argvlen[pos * 2 + 2] = strlen(keys + offset);
            offset += argvlen[pos * 2 + 2];
            argv[pos * 2 + 3] = (char *)(buffer + ((i - *its) * n + v));
            argvlen[pos * 2 + 3] = sizeof(double);
            pos++;
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
}

// hash 表获取一维double物理场
int redis_hmgetd1d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, double *buffer)
{
    int i, j, v;
	int errno = 0;
    int n = *nn;
    int kv_cnt = ((*ite - *its) / *istride + 1) * n;

    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    /*string split*/
    char var[n][32];
    char* token = strtok(varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for (i = *its; i <= *ite; i += *istride) {
        for (v = 0; v < n; v++) {
            sprintf(keys + offset, "%d:%s", i, var[v]);
            argv[pos + 2] = keys + offset;
            argvlen[pos + 2] = strlen(keys + offset);
            offset += argvlen[pos + 2];
            pos++;
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnd(cc, reply, buffer, &mm);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
	return errno;
} 


// hash 表设置二维double物理场
int redis_hmsetd2d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, double *buffer)
{
    int i, j, v;
	int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v =0; v<n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(double);
                pos++;
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
			errno = -1;
			goto err;
    	}
    }

    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
			errno = -1;
			freeReplyObject(reply);
			goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
		if (errno) 
			goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
	return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(double);
                pos++;
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

// hash 表获取二维double物理场
int redis_hmgetd2d(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];
  
    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
    	}
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v=0; v<n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnd(cc, reply, buffer, &mm);
    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

// hash 表设置三维double物理场
int redis_hmsetd(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(double);
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
    	}
    }

    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) 
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(double);
                    pos++;
                }
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

// hash 表获取三维double物理场
int redis_hmgetd(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
    	}
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnd(cc, reply, buffer, &mm);
            
    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

// hash 表完成da二维三维float物理场输出
int redis_da_outputf(redisClusterContext *cc, char *hkey, char *varlist, int *lonids, int *lonide, int *lonstep, 
                      int *latids, int *latide, int *latstep, int *mpas_num_lev_start, int *mpas_num_lev, int *zstep, 
                      int *num_2d, int *num_3d, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j * n;

    int i_cnt = (*lonide - *lonids) / *lonstep + 1;
    int j_cnt = (*latide - *latids) / *latstep + 1;
    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = j_cnt * n;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char*));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *lonids; i <= *lonide; i += *lonstep) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *latids; j <= *latide; j += *latstep) {
            for(v = 0; v < *num_2d; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[2 * pos + 2] = keys + offset;
                argvlen[2 * pos + 2] = strlen(keys + offset);
                offset += argvlen[2 * pos + 2];
                argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                    + (j - *latids) * n + v));
                argvlen[2 * pos + 3] = sizeof(float);
                pos++;
            }
            for(v = 0; v < *num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[*num_2d + v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                        +(j - *latids) * n
                                        +(*num_2d + v * (*mpas_num_lev)) + k - 1));
                    argvlen[2 * pos + 3] = sizeof(float);
                    pos++;
                }
            }
        }
        errno = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
        }
    }

    for(i = *lonids; i <= *lonide; i += *lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;

    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    if (n * sec_i * sec_j < 518400) {
        int kv_cnt = i_cnt * j_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));

        argv[0] = "hmset";
        argvlen[0] = 5;
        argv[1] = hkey;
        argvlen[1] = strlen(hkey);

        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(i = *lonids; i <= *lonide; i += *lonstep) {
            for(j = *latids; j <= *latide; j += *latstep) {
                for(v = 0; v < *num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                        + (j - *latids) * n + v));
                    argvlen[2 * pos + 3] = sizeof(float);
                    pos++;
                }
                for(v = 0; v < *num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(float);
                        pos++;
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
        errno = reply_check(cc, reply);

        free(argv);
        free(argvlen);
        free(keys);
        freeReplyObject(reply);
        return errno;
    } else {
        redisReply *reply;
        int kv_cnt = i_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));
        
        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(j = *latids; j <= *latide; j += *latstep) {
            argv[0] = "hmset";
            argvlen[0] = 5;
            argv[1] = hkey;
            argvlen[1] = strlen(hkey);

            pos = 0;
            offset = 0;
            for(i = *lonids; i <= *lonide; i += *lonstep) {
                for(v=0; v<*num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                        + (j - *latids) * n + v));
                    argvlen[2 * pos + 3] = sizeof(float);
                    pos++;
                }
                for(v = 0; v < *num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(float);
                        pos++;
                    }
                }
            }
            errno = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
            if (errno) {
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                goto err;
            }
        }

        for(j = *latids; j <= *latide; j+=*latstep) {
            int r = redisClusterGetReply(cc, (void **) &reply);
            if (r == REDIS_ERR) {
                printf("Generic Redis Reply Error\n"); 
                errno = -1;
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            if (errno)
                goto err;
        }
    err:
        free(argv);
        free(argvlen);
        free(keys);
        return errno;
    }
#endif
}

// hash 表完成da二维三维double物理场输出
int redis_da_outputd(redisClusterContext *cc, char *hkey, char *varlist, int *lonids, int *lonide, int *lonstep, 
                      int *latids, int *latide, int *latstep, int *mpas_num_lev_start, int *mpas_num_lev, int *zstep, 
                      int *num_2d, int *num_3d, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j * n;

    int i_cnt = (*lonide - *lonids) / *lonstep + 1;
    int j_cnt = (*latide - *latids) / *latstep + 1;
    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = j_cnt * n;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char*));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *lonids; i <= *lonide; i += *lonstep) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *latids; j <= *latide; j += *latstep) {
            for(v = 0; v < *num_2d; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[2 * pos + 2] = keys + offset;
                argvlen[2 * pos + 2] = strlen(keys + offset);
                offset += argvlen[2 * pos + 2];
                argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                    + (j - *latids) * n + v));
                argvlen[2 * pos + 3] = sizeof(double);
                pos++;
            }
            for(v = 0; v < *num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[*num_2d + v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                        + (j - *latids) * n
                                        + (*num_2d + v * (*mpas_num_lev)) + k - 1));
                    argvlen[2 * pos + 3] = sizeof(double);
                    pos++;
                }
            }
        }
        errno = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
        }
    }

    for(i = *lonids; i <= *lonide; i += *lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }
    
err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    if (n * sec_i * sec_j < 518400) {
        int kv_cnt = i_cnt * j_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));

        argv[0] = "hmset";
        argvlen[0] = 5;
        argv[1] = hkey;
        argvlen[1] = strlen(hkey);

        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(i = *lonids; i <= *lonide; i += *lonstep) {
            for(j = *latids; j <= *latide; j += *latstep) {
                for(v = 0; v < *num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                        + (j - *latids) * n + v));
                    argvlen[2 * pos + 3] = sizeof(double);
                    pos++;
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(double);
                        pos++;
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
        errno = reply_check(cc, reply);

        free(argv);
        free(argvlen);
        free(keys);
        freeReplyObject(reply);
        return errno;
    } else {
        redisReply *reply;
        int kv_cnt = i_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));
        
        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(j = *latids; j <= *latide; j += *latstep) {
            argv[0] = "hmset";
            argvlen[0] = 5;
            argv[1] = hkey;
            argvlen[1] = strlen(hkey);

            pos = 0;
            offset = 0;
            for(i = *lonids; i <= *lonide; i += *lonstep) {
                for(v = 0; v < *num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                        + (j - *latids) * n 
                                        + v));
                    argvlen[2 * pos + 3] = sizeof(double);
                    pos++;
                }
                for(v = 0; v < *num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(double);
                        pos++;
                    }
                }
            }
            errno = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
            if (errno) {
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                goto err;
            }
        }

        for(j = *latids; j <= *latide; j += *latstep) {
            int r = redisClusterGetReply(cc, (void **) &reply);
            if (r == REDIS_ERR) { 
                printf("Generic Redis Reply Error\n"); 
                errno = -1;
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            goto err;
        }

    err:
        free(argv);
        free(argvlen);
        free(keys);
        return errno;
    }
#endif
}


// 非阻塞接口
//-------------------------------------------------
// non-block interfaces
//-------------------------------------------------



int redis_hmsetf1d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
    int kv_cnt = ((*ite - *its) / *istride + 1) * n;

    int argc = kv_cnt * 2 + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    /*string split*/
    char var[n][32];
    char* token = strtok(varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for (i = *its; i <= *ite; i += *istride) {
        for (v = 0; v < n; v++) {
            sprintf(keys + offset, "%d:%s", i, var[v]);
            argv[pos * 2 + 2] = keys + offset;
            argvlen[pos * 2 + 2] = strlen(keys + offset);
            offset += argvlen[pos * 2 + 2];
            argv[pos * 2 + 3] = (char *)(buffer + ((i - *its) * n + v));
            argvlen[pos * 2 + 3] = sizeof(float);
            pos++;
        }
    }

    int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    if (ret == REDIS_ERR) {
        printf("append command error\n");
        err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
        errno = -1;
    }

    free(argv);
    free(argvlen);
    free(keys);
    return errno;
}

int redis_hmsetf1d_nonblock_get_reply(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, float *buffer) {
    void *reply;
    int errno = 0;
    int status = redisClusterGetReply(cc, &reply);
    if (status == REDIS_OK) {
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
    }
    else {
        printf("redis_hmsetf1d_nonblock_get_reply error!\n");
        err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
        errno = -1;
    }
    return errno;
}

int redis_hmgetf1d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
    int kv_cnt = ((*ite - *its) / *istride + 1) * n;

    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    /*string split*/
    char var[n][32];
    char* token = strtok(varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for (i = *its; i <= *ite; i += *istride) {
        for (v = 0; v < n; v++) {
            sprintf(keys + offset, "%d:%s", i, var[v]);
            argv[pos + 2] = keys + offset;
            argvlen[pos + 2] = strlen(keys + offset);
            offset += argvlen[pos + 2];
            pos++;
        }
    }

    int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    if (ret == REDIS_ERR) {
        printf("append command error\n");
        err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
        errno = -1;
    }

    free(argv);
    free(argvlen);
    free(keys);
    return errno;
} 

int redis_hmgetf1d_nonblock_get_reply(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *nn, float *buffer) {
    void *reply;
    int errno = 0;
    int status = redisClusterGetReply(cc, &reply);
    if (status == REDIS_OK) {
        int mm = 0;        
        errno = reply_returnf(cc, reply, buffer, &mm);
        err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
        freeReplyObject(reply);
    }
    else {
        printf("redis_hmsetf1d_nonblock_get_reply error!\n");
        err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
        errno = -1;
    }
    return errno;
}

int redis_hmsetf2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v =0; v<n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(float);
                pos++;
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }
    

    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }
    
err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(float);
                pos++;
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_hmgetf2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];
  
    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
	
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
		int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v=0; v<n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnf(cc, reply, buffer, &mm);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_hmsetf_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(float);
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(float);
                    pos++;
                }
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_hmgetf_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
		offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
 		int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnf(cc, reply, buffer, &mm);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_hmsetd2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v =0; v<n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(double);
                pos++;
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos * 2 + 2] = keys + offset;
                argvlen[pos * 2 + 2] = strlen(keys + offset);
                offset += argvlen[pos * 2 + 2];
                argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jn
                                            + (j - *jts) * n + v));
                argvlen[pos * 2 + 3] = sizeof(double);
                pos++;
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_hmgetd2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int jn_cnt = j_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];
  
    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v = 0; v < n; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }
	
	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(v=0; v<n; v++) {
                sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                argv[pos + 2] = keys + offset;
                argvlen[pos + 2] = strlen(keys + offset);
                offset += argvlen[pos + 2];
                pos++;
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnd(cc, reply, buffer, &mm);
    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_hmsetd_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(double);
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmset";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos * 2 + 2] = keys + offset;
                    argvlen[pos * 2 + 2] = strlen(keys + offset);
                    offset += argvlen[pos * 2 + 2];
                    argv[pos * 2 + 3] = (char *)(buffer+((i - *its) * sec_jkn
                                                + (j - *jts) * sec_kn 
                                                + (k - *kts) * n
                                                + v));
                    argvlen[pos * 2 + 3] = sizeof(double);
                    pos++;
                }
            }
        }
    }

    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_check(cc, reply);

    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_hmgetd_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, int *ite, int *istride,
                  int *jts, int *jte, int *jstride, int *kts, int *kte, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kte - *kts + 1;
    int sec_jkn = sec_j * sec_k * n;
    int sec_kn = sec_k * n;
    
    int i_cnt = (*ite - *its) / *istride + 1;
    int j_cnt = (*jte - *jts) / *jstride + 1;
    int k_cnt = (*kte - *kts) / *kstride + 1;
    int jkn_cnt = j_cnt * k_cnt * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));

    for(i = *its; i <= *ite; i += *istride) {
        argv[0] = "hmget";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *jts; j <= *jte; j += *jstride) {
            for (k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i += *istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    int kv_cnt = i_cnt * jkn_cnt;
    int argc = kv_cnt + 2;
    char** argv = calloc(kv_cnt + 2, sizeof(char *));
    size_t* argvlen = calloc(kv_cnt + 2, sizeof(size_t));

    argv[0] = "hmget";
    argvlen[0] = 5;
    argv[1] = hkey;
    argvlen[1] = strlen(hkey);

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *its; i <= *ite; i += *istride) {
        for(j = *jts; j <= *jte; j += *jstride) {
            for(k = *kts; k <= *kte; k += *kstride) {
                for(v = 0; v < n; v++) {
                    sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[v]);
                    argv[pos + 2] = keys + offset;
                    argvlen[pos + 2] = strlen(keys + offset);
                    offset += argvlen[pos + 2];
                    pos++;
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
    errno = reply_returnd(cc, reply, buffer, &mm);
            
    free(argv);
    free(argvlen);
    free(keys);
    freeReplyObject(reply);
    return errno;
#endif
}

int redis_da_outputf_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *lonids, int *lonide, int *lonstep, 
                      int *latids, int *latide, int *latstep, int *mpas_num_lev_start, int *mpas_num_lev, int *zstep, 
                      int *num_2d, int *num_3d, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j * n;

    int i_cnt = (*lonide - *lonids) / *lonstep + 1;
    int j_cnt = (*latide - *latids) / *latstep + 1;
    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = j_cnt * n;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char*));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *lonids; i <= *lonide; i += *lonstep) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *latids; j <= *latide; j += *latstep) {
            for(v = 0; v < *num_2d; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[2 * pos + 2] = keys + offset;
                argvlen[2 * pos + 2] = strlen(keys + offset);
                offset += argvlen[2 * pos + 2];
                argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                    + (j - *latids) * n + v));
                argvlen[2 * pos + 3] = sizeof(float);
                pos++;
            }
            for(v = 0; v < *num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[*num_2d + v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                        +(j - *latids) * n
                                        +(*num_2d + v * (*mpas_num_lev)) + k - 1));
                    argvlen[2 * pos + 3] = sizeof(float);
                    pos++;
                }
            }
        }
        int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
    	if (ret == REDIS_ERR) {
       		printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
    	}
    }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *lonids; i <= *lonide; i += *lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    if (n * sec_i * sec_j < 518400) {
        int kv_cnt = i_cnt * j_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));

        argv[0] = "hmset";
        argvlen[0] = 5;
        argv[1] = hkey;
        argvlen[1] = strlen(hkey);

        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(i = *lonids; i <= *lonide; i += *lonstep) {
            for(j = *latids; j <= *latide; j += *latstep) {
                for(v = 0; v < *num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                        + (j - *latids) * n + v));
                    argvlen[2 * pos + 3] = sizeof(float);
                    pos++;
                }
                for(v = 0; v < *num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(float);
                        pos++;
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
        errno = reply_check(cc, reply);

        free(argv);
        free(argvlen);
        free(keys);
        freeReplyObject(reply);
        return errno;
    } else {
        redisReply *reply;
        int kv_cnt = i_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));
        
        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(j = *latids; j <= *latide; j += *latstep) {
            argv[0] = "hmset";
            argvlen[0] = 5;
            argv[1] = hkey;
            argvlen[1] = strlen(hkey);

            pos = 0;
            offset = 0;
            for(i = *lonids; i <= *lonide; i += *lonstep) {
                for(v=0; v<*num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                        + (j - *latids) * n + v));
                    argvlen[2 * pos + 3] = sizeof(float);
                    pos++;
                }
                for(v = 0; v < *num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(float);
                        pos++;
                    }
                }
            }
            int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
			if (ret == REDIS_ERR) {
				printf("append command error\n");
                errno = -1;
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                goto err;
			}
        }

		errno = redis_sendall(cc);
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
            goto err;
        }

        for(j = *latids; j <= *latide; j+=*latstep) {
            int r = redisClusterGetReply(cc, (void **) &reply);
            if (r == REDIS_ERR) { 
                printf("Generic Redis Reply Error\n"); 
                errno = -1;
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            if (errno)
                goto err;
        }

err:
        free(argv);
        free(argvlen);
        free(keys);
        return errno;
    }
#endif
}

int redis_da_outputd_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *lonids, int *lonide, int *lonstep, 
                      int *latids, int *latide, int *latstep, int *mpas_num_lev_start, int *mpas_num_lev, int *zstep, 
                      int *num_2d, int *num_3d, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j * n;

    int i_cnt = (*lonide - *lonids) / *lonstep + 1;
    int j_cnt = (*latide - *latids) / *latstep + 1;
    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
    redisReply *reply;
    int kv_cnt = j_cnt * n;
    int argc = 2 * kv_cnt + 2;
    char** argv = calloc(2 * kv_cnt + 2, sizeof(char*));
    size_t* argvlen = calloc(2 * kv_cnt + 2, sizeof(size_t));
    char lon_hkey[MAX_KEY_LEN];

    int pos = 0, offset = 0;
    char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
    for(i = *lonids; i <= *lonide; i += *lonstep) {
        argv[0] = "hmset";
        argvlen[0] = 5;
#ifdef LON_MAP
        sprintf(lon_hkey, "%s{%d}", hkey, i);
#else
        sprintf(lon_hkey, "%s:%d", hkey, i);
#endif
        argv[1] = lon_hkey;
        argvlen[1] = strlen(lon_hkey);

        pos = 0;
        offset = 0;
        for(j = *latids; j <= *latide; j += *latstep) {
            for(v = 0; v < *num_2d; v++) {
                sprintf(keys + offset, "%d:0:%s", j, var[v]);
                argv[2 * pos + 2] = keys + offset;
                argvlen[2 * pos + 2] = strlen(keys + offset);
                offset += argvlen[2 * pos + 2];
                argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                    + (j - *latids) * n + v));
                argvlen[2 * pos + 3] = sizeof(double);
                pos++;
            }
            for(v = 0; v < *num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                    sprintf(keys + offset, "%d:%d:%s", j, k, var[*num_2d + v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                        + (j - *latids) * n
                                        + (*num_2d + v * (*mpas_num_lev)) + k - 1));
                    argvlen[2 * pos + 3] = sizeof(double);
                    pos++;
                }
            }
        }
		int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
		if (ret == REDIS_ERR) {
			printf("append command error\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            goto err;
		}
    }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *lonids; i <= *lonide; i += *lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            errno = -1;
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno)
            goto err;
    }

err:
    free(argv);
    free(argvlen);
    free(keys);
    return errno;
#else
    if (n * sec_i * sec_j < 518400) {
        int kv_cnt = i_cnt * j_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));

        argv[0] = "hmset";
        argvlen[0] = 5;
        argv[1] = hkey;
        argvlen[1] = strlen(hkey);

        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(i = *lonids; i <= *lonide; i += *lonstep) {
            for(j = *latids; j <= *latide; j += *latstep) {
                for(v = 0; v < *num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                        + (j - *latids) * n + v));
                    argvlen[2 * pos + 3] = sizeof(double);
                    pos++;
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer+((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(double);
                        pos++;
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommandArgv(cc, argc, (const char **)argv, argvlen);
        errno = reply_check(cc, reply);

        free(argv);
        free(argvlen);
        free(keys);
        freeReplyObject(reply);
        return errno;
    } else {
        redisReply *reply;
        int kv_cnt = i_cnt * n;
        int argc = 2 * kv_cnt + 2;
        char** argv = calloc(argc, sizeof(char *));
        size_t* argvlen = calloc(argc, sizeof(size_t));
        
        int pos = 0, offset = 0;
        char* keys = calloc(kv_cnt * MAX_FIELD_LEN, sizeof(char));
        for(j = *latids; j <= *latide; j += *latstep) {
            argv[0] = "hmset";
            argvlen[0] = 5;
            argv[1] = hkey;
            argvlen[1] = strlen(hkey);

            pos = 0;
            offset = 0;
            for(i = *lonids; i <= *lonide; i += *lonstep) {
                for(v = 0; v < *num_2d; v++) {
                    sprintf(keys + offset, "%d:%d:0:%s", i, j, var[v]);
                    argv[2 * pos + 2] = keys + offset;
                    argvlen[2 * pos + 2] = strlen(keys + offset);
                    offset += argvlen[2 * pos + 2];
                    argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                        + (j - *latids) * n 
                                        + v));
                    argvlen[2 * pos + 3] = sizeof(double);
                    pos++;
                }
                for(v = 0; v < *num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k += *zstep) {
                        sprintf(keys + offset, "%d:%d:%d:%s", i, j, k, var[*num_2d + v]);
                        argv[2 * pos + 2] = keys + offset;
                        argvlen[2 * pos + 2] = strlen(keys + offset);
                        offset += argvlen[2 * pos + 2];
                        argv[2 * pos + 3] = (char *)(buffer + ((i - *lonids) * sec_kjn
                                            +(j - *latids) * n
                                            +(*num_2d + v*(*mpas_num_lev)) + k - 1));
                        argvlen[2 * pos + 3] = sizeof(double);
                        pos++;
                    }
                }
            }
            int ret = redisClusterAppendCommandArgv(cc, argc, (const char **)argv, argvlen);
			if (ret == REDIS_ERR) {
				printf("append command error\n");
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                errno = -1;
                goto err;
			}
        }

		errno = redis_sendall(cc);
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
            goto err;
        }

        for(j = *latids; j <= *latide; j += *latstep) {
            int r = redisClusterGetReply(cc, (void **) &reply);
            if (r == REDIS_ERR) { 
                printf("Generic Redis Reply Error\n"); 
                errno = -1;
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            if (errno)
                goto err;
        }
err:
        free(argv);
        free(argvlen);
        free(keys);
        return errno;
    }
#endif
}


//----------------------------------------------------------------
// reply_interface
//----------------------------------------------------------------
// 处理接受值为 double
int reply_returnd(redisClusterContext *cc, redisReply *reply, double *buf, int *tag)
{
    int i;
	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}
    switch(reply->type) {
        case 1 :
            memcpy(buf, reply->str, sizeof(double));
            break;
        case 2:
            for(i = 0; i < reply->elements; i++) {
                memcpy(buf + *tag, reply->element[i]->str, sizeof(double));
                (*tag)++;
            }
            break;
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;
    }
	return errno;
}

// 处理接受值为 float
int reply_returnf(redisClusterContext *cc, redisReply *reply, float *buf, int *tag)
{
    int i;
	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}
    switch(reply->type) {
        case 1 :
            memcpy(buf, reply->str, sizeof(float));
            break;
        case 2:
            for(i = 0; i < reply->elements; i++) {
                memcpy(buf + *tag, reply->element[i]->str, sizeof(float));
                (*tag)++;
            }
            break;
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;    }
	return errno;
}

// 处理接受值为 int
int reply_returni(redisClusterContext *cc, redisReply *reply, int *buf, int *tag)
{
    int i;
	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}
    switch(reply->type) {
        case 1 :
            memcpy(buf, reply->str, sizeof(int));
            break;
        case 2:
            for(i = 0; i < reply->elements; i++) {
                memcpy(buf + *tag, reply->element[i]->str, sizeof(int));
                (*tag)++;
            }
            break;
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;
    }
	return errno;
}

// 处理接受值为 string
int reply_returns(redisClusterContext *cc, redisReply *reply, char *buf, int *len) {
	int i;
	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}
    switch(reply->type) {
        case 1 :
			memcpy(buf, reply->str, reply->len);
			*len = reply->len;
            break;
        case 2:
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;
    }
	return errno;
}

// 不需要接受值 
int reply_check(redisClusterContext *cc, redisReply *reply)
{
    if(NULL == reply) {
        err_print("Get Null Reply:%s", cc->errstr);
		return -1;
    }
	return 0;
}

// 字符串接口
#else
// str interface

#include "sb.h"

int RedisHgeti(redisClusterContext *cc, char *hash, char *key, int *value)
{
    redisReply *reply = redisClusterCommand(cc, "HGET %s %s", hash, key);
    int errno = reply_returni(cc, reply, value, 0);
    freeReplyObject(reply);
    return errno;
}

int RedisHgetf(redisClusterContext *cc, char *hash, char *key, float *value)
{
    redisReply *reply = redisClusterCommand(cc, "HGET %s %s", hash, key);
    int errno = reply_returnf(cc, reply, value, 0);
    freeReplyObject(reply);
    return errno;
}

int RedisHgetd(redisClusterContext *cc, char *hash, char *key, double *value)
{
    redisReply *reply = redisClusterCommand(cc, "HGET %s %s", hash, key);
    int errno = reply_returnd(cc, reply, value, 0);
    freeReplyObject(reply);
    return errno;
}

int RedisHgets(redisClusterContext *cc, char *hash, char *key, 
        char *value, int *len1, int *len2, int *len3)
{
    redisReply *reply = redisClusterCommand(cc, "HGET %s %s", hash, key);
    int errno = reply_returns(cc, reply, value, 0);
    freeReplyObject(reply);
    return errno;
}


int RedisHseti(redisClusterContext *cc, char *hash, char *key, int *value)
{
    redisReply *reply = redisClusterCommand(cc, "HSET %s %s %d", hash, key, *value);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
    return errno;
}

int RedisHsetf(redisClusterContext *cc, char *hash, char *key, float *value)
{
    redisReply *reply = redisClusterCommand(cc, "HSET %s %s %f", hash, key, *value);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
    return errno;
}

int RedisHsetd(redisClusterContext *cc, char *hash, char *key, double *value)
{
    redisReply *reply = redisClusterCommand(cc, "HSET %s %s %.12f", hash, key, *value);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
    return errno;
}

int RedisHsets(redisClusterContext *cc, char *hash, char *key, char *value)
{
    redisReply *reply = redisClusterCommand(cc, "HSET %s %s %s", hash, key, value);
    int errno = reply_check(cc, reply);
    freeReplyObject(reply);
    return errno;
}

int cmpfunc(const void * a, const void * b)
{
   return (*(int*)a - *(int*)b);
}

int cmpf(const void *a, const void *b)
{
    struct obsF fa, fb;
    fa = *(struct obsF*)a;
    fb = *(struct obsF*)b;
    float dx,dy,dlev;
    dx = (fa.lon) - (fb.lon);
    if (dx < 1e-6 && dx > -1e-6) {
        dy = (fa.lat) - (fb.lat);
        if (dy < 1e-6 && dy > -1e-6) {
            dlev = (fa.lev) - (fb.lev);
            if (dlev > 0) return 1;
            else return -1;
        }
        else {
            if (dy > 0) return 1;
            else return -1;
        }
    }
    else {
        if (dx > 0) return 1;
        else return -1;
    }
}

int cmpd(const void *a, const void *b)
{
    struct obsD fa, fb;
    fa = *(struct obsD*)a;
    fb = *(struct obsD*)b;
    double dx,dy,dlev;
    dx = (fa.lon) - (fb.lon);
    if (dx < 1e-6 && dx > -1e-6) {
        dy = (fa.lat) - (fb.lat);
        if (dy < 1e-6 && dy > -1e-6) {
            dlev = (fa.lev) - (fb.lev);
            if (dlev > 0) return 1;
            else return -1;
        }
        else {
            if (dy > 0) return 1;
            else return -1;
        }
    }
    else {
        if (dx > 0) return 1;
        else return -1;
    }
}


int redis_obscount(redisClusterContext *cc, char *hkey, float *its, float *ite, float *jts, float *jte, int *buf)
{
    int errno = 0;
    redisReply *reply;
    StringBuilder *key = sb_create();
    int i;

    //range lon
    sb_appendf(key, "zrangebyscore lon_index:%s %.5f %.5f ", hkey, (*its), (*ite));
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_ARRAY) {
        printf("Wrong zrange type\n");
        return -1;
    }
    int num_lon;
    int *lon_indx;
    num_lon = reply->elements;
    lon_indx = (int*)malloc(sizeof(int)*num_lon);
    for (i = 0; i < num_lon; i++) {
        lon_indx[i] = atoi(reply->element[i]->str);
    }
    freeReplyObject(reply);
    sb_reset(key);

    //range lat
    sb_appendf(key, "zrangebyscore lat_index:%s %.5f %.5f ", hkey, (*jts), (*jte));
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_ARRAY) {
        printf("Wrong zrange type\n");
        return -1;
    }
    int num_lat;
    int *lat_indx;
    num_lat = reply->elements;
    lat_indx = (int*)malloc(sizeof(int)*num_lat);
    for (i = 0; i < num_lat; i++) {
        lat_indx[i] = atoi(reply->element[i]->str);
    }
    freeReplyObject(reply);
    sb_reset(key);

    int *indx;
    indx = (int*)malloc(sizeof(int)*(num_lon+num_lat));
    for (i = 0; i < num_lon; i++) indx[i] = lon_indx[i];
    for (i = 0; i < num_lat; i++) indx[num_lon+i] = lat_indx[i];
    qsort(indx, num_lat+num_lon, sizeof(int), cmpfunc);

    int cnt;
    cnt = 0;
    for (i = 1; i < num_lat+num_lon; i++) {
        if (indx[i] == indx[i-1]) {
            cnt++;
        }
    }

    (*buf) = cnt;

    free(lon_indx);
    free(lat_indx);
    free(indx);
    sb_free(key);
    return errno;

}

int redis_getobservef(redisClusterContext *cc, char *hkey, char *varlist, float *its, float *ite, float *jts, float *jte,
                        float *lon_buf, float *lat_buf, float *lev_buf, float *val_buf, int *nn)
{
    int errno = 0;
    redisReply *reply;
    StringBuilder *key = sb_create();
    int i;

    //range lon
    sb_appendf(key, "zrangebyscore lon_index:%s %.5f %.5f ", hkey, (*its), (*ite));
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_ARRAY) {
        printf("Wrong zrange type\n");
        return -1;
    }
    int num_lon;
    int *lon_indx;
    num_lon = reply->elements;
    lon_indx = (int*)malloc(sizeof(int)*num_lon);
    for (i = 0; i < num_lon; i++) {
        lon_indx[i] = atoi(reply->element[i]->str);
    }
    freeReplyObject(reply);
    sb_reset(key);

    //range lat
    sb_appendf(key, "zrangebyscore lat_index:%s %.5f %.5f ", hkey, (*jts), (*jte));
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_ARRAY) {
        printf("Wrong zrange type\n");
        return -1;
    }
    int num_lat;
    int *lat_indx;
    num_lat = reply->elements;
    lat_indx = (int*)malloc(sizeof(int)*num_lat);
    for (i = 0; i < num_lat; i++) {
        lat_indx[i] = atoi(reply->element[i]->str);
    }
    freeReplyObject(reply);
    sb_reset(key);

    int *indx, *indx_r;
    indx = (int*)malloc(sizeof(int)*(num_lon+num_lat));
    indx_r = (int*)malloc(sizeof(int)*(num_lon+num_lat));
    for (i = 0; i < num_lon; i++) indx[i] = lon_indx[i];
    for (i = 0; i < num_lat; i++) indx[num_lon+i] = lat_indx[i];
    qsort(indx, num_lat+num_lon, sizeof(int), cmpfunc);

    int cnt;
    cnt = 0;
    for (i = 1; i < num_lat+num_lon; i++) {
        if (indx[i] == indx[i-1]) {
            indx_r[cnt] = indx[i];
            cnt++;
        }
    }

    //printf("cnt:%d\n", cnt);

    //string split
    int n = (*nn);

    char var[n][32];
    char* token = strtok( varlist, ":");

    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    int j;
    int mm;

    int lim; //batch_insert
    lim = 100000;
    int st, ed;
    st = 0;

    //hmget lon
    if (cnt > 0) {
        for (st = 0; st < cnt; st+=lim) {
            ed = st + lim;
            if (ed > cnt) ed = cnt;

            sb_appendf(key, "hmget lon:%s ", hkey);
            for (i = st; i < ed; i++) {
                sb_appendf(key, "%d ", indx_r[i]); 
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnf(cc, reply, lon_buf + st, &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);

            if (cnt == 216915) {
            //    printf("get_lon\n");
            }

            //hmget lat
            sb_appendf(key, "hmget lat:%s ", hkey);
            for (i = st; i < ed; i++) {
                sb_appendf(key, "%d ", indx_r[i]); 
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnf(cc, reply, lat_buf + st, &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);

            if (cnt == 216915) {
            //    printf("get_lat\n");
            }

            //hmget lev
            sb_appendf(key, "hmget lev:%s ", hkey);
            for (i = st; i < ed; i++) {
                sb_appendf(key, "%d ", indx_r[i]); 
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnf(cc, reply, lev_buf + st, &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);

            if (cnt == 216915) {
            //    printf("get_lev\n");
            }

            //hmget lon
            sb_appendf(key, "hmget val:%s ", hkey);
            for (i = st; i < ed; i++) {
                for (j = 0; j < n; j++) {
                    sb_appendf(key, "%d:%s ", indx_r[i], var[j]); 
                }
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnf(cc, reply, val_buf + (st * n), &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);
            if (cnt == 216915) {
             //   printf("get_val\n");
            }

        }
    }

    free(lon_indx);
    free(lat_indx);
    free(indx);
    free(indx_r);
    sb_free(key);
    return errno;

}

int redis_getobserved(redisClusterContext *cc, char *hkey, char *varlist, double *its, double *ite, double *jts, double *jte,
                        double *lon_buf, double *lat_buf, double *lev_buf, double *val_buf, int *nn)
{
    int errno = 0;
    redisReply *reply;
    StringBuilder *key = sb_create();
    int i;

    //range lon
    sb_appendf(key, "zrangebyscore lon_index:%s %.6lf %.6lf ", hkey, (*its), (*ite));
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_ARRAY) {
        printf("Wrong zrange type\n");
        return -1;
    }
    int num_lon;
    int *lon_indx;
    num_lon = reply->elements;
    lon_indx = (int*)malloc(sizeof(int)*num_lon);
    for (i = 0; i < num_lon; i++) {
        lon_indx[i] = atoi(reply->element[i]->str);
    }
    freeReplyObject(reply);
    sb_reset(key);

    //range lat
    sb_appendf(key, "zrangebyscore lat_index:%s %.6lf %.6lf ", hkey, (*jts), (*jte));
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_ARRAY) {
        printf("Wrong zrange type\n");
        return -1;
    }
    int num_lat;
    int *lat_indx;
    num_lat = reply->elements;
    lat_indx = (int*)malloc(sizeof(int)*num_lat);
    for (i = 0; i < num_lat; i++) {
        lat_indx[i] = atoi(reply->element[i]->str);
    }
    freeReplyObject(reply);
    sb_reset(key);

    int *indx, *indx_r;
    indx = (int*)malloc(sizeof(int)*(num_lon+num_lat));
    indx_r = (int*)malloc(sizeof(int)*(num_lon+num_lat));
    for (i = 0; i < num_lon; i++) indx[i] = lon_indx[i];
    for (i = 0; i < num_lat; i++) indx[num_lon+i] = lat_indx[i];
    qsort(indx, num_lat+num_lon, sizeof(int), cmpfunc);

    int cnt;
    cnt = 0;
    for (i = 1; i < num_lat+num_lon; i++) {
        if (indx[i] == indx[i-1]) {
            indx_r[cnt] = indx[i];
            cnt++;
        }
    }

    //printf("cnt:%d\n", cnt);

    //string split
    int n = (*nn);

    char var[n][32];
    char* token = strtok( varlist, ":");

    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    int j;
    int mm;
    //hmget lon
    int lim; //batch_insert
    lim = 100000;
    int st, ed;
    st = 0;

    if (cnt > 0) {
        for (st = 0; st < cnt; st+=lim) {
            ed = st + lim;
            if (ed > cnt) ed = cnt;
            sb_appendf(key, "hmget lon:%s ", hkey);
            for (i = st; i < ed; i++) {
                sb_appendf(key, "%d ", indx_r[i]); 
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnd(cc, reply, lon_buf + st, &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);

            //hmget lat
            sb_appendf(key, "hmget lat:%s ", hkey);
            for (i = st; i < ed; i++) {
                sb_appendf(key, "%d ", indx_r[i]); 
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnd(cc, reply, lat_buf + st, &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);

            //hmget lev
            sb_appendf(key, "hmget lev:%s ", hkey);
            for (i = st; i < ed; i++) {
                sb_appendf(key, "%d ", indx_r[i]); 
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnd(cc, reply, lev_buf + st, &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);

            //hmget lon
            sb_appendf(key, "hmget val:%s ", hkey);
            for (i = st; i < ed; i++) {
                for (j = 0; j < n; j++) {
                    sb_appendf(key, "%d:%s ", indx_r[i], var[j]); 
                }
            }
            mm = 0;
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_returnd(cc, reply, val_buf + (st * n), &mm);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);
        }
    }

    free(lon_indx);
    free(lat_indx);
    free(indx);
    free(indx_r);
    sb_free(key);
    return errno;

}
// 非规则观测接口，输入观测
// zset (lon/lat_index):hkey
// hset (lon/lat/lev/val:hkey)
int redis_setobservef(redisClusterContext *cc, char *hkey, char *varlist, float *lon_buf, float *lat_buf, 
                        float *lev_buf, float *val_buf, int *obs_num, int *nn)
{
    int errno = 0;
    redisReply *reply;
    StringBuilder *key = sb_create();
    int i;

    //string split
    int n = (*nn);
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    // zcard
    sb_appendf(key, "zcard lon_index:%s ", hkey);
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_INTEGER) return -1;
    int num_obs;
    num_obs = (int)(reply->integer);
    freeReplyObject(reply);
    sb_reset(key);

    //sort
    struct obsF *obs;
    int *p;
    float dx,dy,dlev;
    obs = (struct obsF*)malloc(sizeof(struct obsF)*(*obs_num));
    p = (int*)malloc(sizeof(int)*(*obs_num));
    for (i = 0; i < (*obs_num); i++) {
        obs[i].lon = lon_buf[i];
        obs[i].lat = lat_buf[i];
        obs[i].lev = lev_buf[i];
        obs[i].id = i;
        p[i] = 1;
    }
    qsort(obs, (*obs_num), sizeof(struct obsF), cmpf);
    for (i = 0; i < (*obs_num) - 1; i++) {
        dy = obs[i+1].lat - obs[i].lat;
        dx = obs[i+1].lon - obs[i].lon;
        dlev = obs[i+1].lev - obs[i].lev;
        if (dx < 1e-6 && dx > -1e-6 && dy < 1e-6 && dy > -1e-6 && dlev < 1e-6 && dlev > -1e-6) {
            p[obs[i+1].id] = 0;
        }
    }

    //for_debug!!!
    //num_obs = 0;

    //zadd lon_index
    int lim; //batch_insert
    lim = 100000;
    int st, ed;
    st = 0;

    for (st = 0; st < (*obs_num); st+=lim) {

        ed = st + lim;
        if (ed > (*obs_num)) {
            ed = (*obs_num);
        }

        sb_reset(key);
        sb_appendf(key, "zadd lon_index:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%.6f %d ", lon_buf[i], num_obs+i);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //zadd lat_index
        sb_appendf(key, "zadd lat_index:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%.6f %d ", lat_buf[i], num_obs+i);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset lon
        sb_appendf(key, "hmset lon:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%d %.6f ", num_obs+i, lon_buf[i]);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset lat
        sb_appendf(key, "hmset lat:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%d %.6f ", num_obs+i, lat_buf[i]);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset lev
        sb_appendf(key, "hmset lev:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%d %.6f ", num_obs+i, lev_buf[i]);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset val
        int j;
        for (j = 0; j < n; j++) {
            sb_appendf(key, "hmset val:%s ", hkey);
            for (i = st; i < ed; i++) {
                if (p[i] == 1) {
                    sb_appendf(key, "%d:%s %.6f ", num_obs+i, var[j], val_buf[i*n+j]);
                }
            }
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_check(cc, reply);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);
        }
    }

    sb_free(key);
    free(p);
    free(obs);
    return errno;
}

int redis_setobserved(redisClusterContext *cc, char *hkey, char *varlist, double *lon_buf, double *lat_buf, 
                        double *lev_buf, double *val_buf, int *obs_num, int *nn)
{
    int errno = 0;
    redisReply *reply;
    StringBuilder *key = sb_create();
    int i;

    //string split
    int n = (*nn);
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    // zcard
    sb_appendf(key, "zcard lon_index:%s ", hkey);
    reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    if (errno != 0) return errno;
    if (reply->type != REDIS_REPLY_INTEGER) return -1;
    int num_obs;
    num_obs = (int)reply->integer;
    freeReplyObject(reply);
    sb_reset(key);

    //sort
    struct obsD *obs;
    int *p;
    double dx,dy,dlev;
    obs = (struct obsD*)malloc(sizeof(struct obsD)*(*obs_num));
    p = (int*)malloc(sizeof(int)*(*obs_num));
    for (i = 0; i < (*obs_num); i++) {
        obs[i].lon = lon_buf[i];
        obs[i].lat = lat_buf[i];
        obs[i].lev = lev_buf[i];
        obs[i].id = i;
        p[i] = 1;
    }
    qsort(obs, (*obs_num), sizeof(struct obsD), cmpd);
    for (i = 0; i < (*obs_num) - 1; i++) {
        dy = obs[i+1].lat - obs[i].lat;
        dx = obs[i+1].lon - obs[i].lon;
        dlev = obs[i+1].lev - obs[i].lev;
        if (dx < 1e-6 && dx > -1e-6 && dy < 1e-6 && dy > -1e-6 && dlev < 1e-6 && dlev > -1e-6) {
            p[obs[i+1].id] = 0;
        }
    }

    //for_debug!!!
    //num_obs = 0;
    int lim = 100000;
    int st,ed;

    st = 0;

    for (st = 0; st < (*obs_num); st+=lim) {

        ed = st + lim;
        if (ed > (*obs_num)) {
            ed = (*obs_num);
        }
        sb_reset(key);
        //zadd lon_index
        sb_appendf(key, "zadd lon_index:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%.6lf %d ", lon_buf[i], num_obs+i);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //zadd lat_index
        sb_appendf(key, "zadd lat_index:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%.6lf %d ", lat_buf[i], num_obs+i);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset lon
        sb_appendf(key, "hmset lon:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%d %.6lf ", num_obs+i, lon_buf[i]);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset lat
        sb_appendf(key, "hmset lat:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%d %.6lf ", num_obs+i, lat_buf[i]);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset lev
        sb_appendf(key, "hmset lev:%s ", hkey);
        for (i = st; i < ed; i++) {
            if (p[i] == 1) {
                sb_appendf(key, "%d %.6lf ", num_obs+i, lev_buf[i]);
            }
        }
        reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno != 0) return errno;
        freeReplyObject(reply);
        sb_reset(key);

        //hmset val
        int j;
        for (j = 0; j < n; j++) {
            sb_appendf(key, "hmset val:%s ", hkey);
            for (i = st; i < ed; i++) {
                if (p[i] == 1) {
                    sb_appendf(key, "%d:%s %.6lf ", num_obs+i, var[j], val_buf[i*n+j]);
                }
            }
            reply = redisClusterCommand(cc, sb_concat(key));
            errno = reply_check(cc, reply);
            if (errno != 0) return errno;
            freeReplyObject(reply);
            sb_reset(key);
        }
    }
    
    sb_free(key);
    free(obs);
    free(p);
    return errno;
}

int redis_hmsetf1d(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int size = (*ite - *its + 1)* n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
          for(v=0; v<n; v++) {
              sb_appendf(key, "%d:%s %.6f ", i, var[v], *(buffer+((i - *its) * n + v)));
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
    sb_free(key);
    return errno;
}

int redis_hmsetf2d(redisClusterContext *cc, char *hkey, char *varlist, int *its, 
        int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn + (j - *jts) * n + v)));
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn + (j - *jts) * n + v)));
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:%d:0:%s %.6f ", i, j, var[v], 
                          *(buffer+((i - *its) * sec_jn + (j - *jts) * n + v)));
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
    if (errno) {
        goto err;
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmsetf(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s %.6f ", i, j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
    if (errno) {
        goto err;
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetf1d(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *nn, float *buffer)
{
    int i, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int size = (*ite - *its + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(v=0; v<n; v++) {
            sb_appendf(key, "%d:%s ", i, var[v]); 
        }
    }

    int mm=0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnf(cc, reply, buffer, &mm);
    freeReplyObject(reply);
    sb_free(key);
    return errno;
}

int redis_hmgetf2d(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:%d:0:%s ", i, j, var[v]); 
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnf(cc, reply, buffer, &mm);
    freeReplyObject(reply);
    if (errno) {
        goto err;
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetf(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            exit(-1);
        }
        int errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s ", i, j, k, var[v]); 
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnf(cc, reply, buffer, &mm);
    freeReplyObject(reply);
    if (errno) {
        goto err;
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmsetd1d(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int size = (*ite - *its + 1)* n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
          for(v=0; v<n; v++) {
              sb_appendf(key, "%d:%s %.12f ", i, var[v], *(buffer+((i - *its) * n + v)));
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
    sb_free(key);
    return errno;
}

int redis_hmsetd2d(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn
                          +(j - *jts) * n + v)));
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn
                          +(j - *jts) * n + v)));
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        if (errno) {
            goto err;
        }
        freeReplyObject(reply);
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:%d:0:%s %.12f ", i, j, var[v], 
                          *(buffer+((i - *its) * sec_jn
                          +(j - *jts) * n + v)));
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmsetd(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s %.12f ", i, j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetd1d(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *nn, double *buffer)
{
    int i, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int size = (*ite - *its + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(v=0; v<n; v++) {
            sb_appendf(key, "%d:%s ", i, var[v]); 
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnd(cc, reply, buffer, &mm);
    freeReplyObject(reply);
    sb_free(key);
    return errno;
}

int redis_hmgetd2d(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:%d:0:%s ", i, j, var[v]); 
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnd(cc, reply, buffer, &mm);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetd(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s ", i, j, k, var[v]); 
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnd(cc, reply, buffer, &mm);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_da_outputd(redisClusterContext *cc, char *hkey, char *varlist, 
        int *lonids, int *lonide, int *lonstep, 
        int *latids, int *latide, int *latstep, int *mpas_num_lev_start, 
        int *mpas_num_lev, int *zstep, int *num_2d, int *num_3d, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
	StringBuilder	*key = sb_create();
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j*n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }
#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#endif
#else
    if (n * sec_i * sec_j < 518400) {
        sb_appendf(key, "hmset %s ", hkey);
        for(i = *lonids; i <= *lonide; i+=*lonstep) {
            for(j = *latids; j <= *latide; j+=*latstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.12f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.12f ", i, j, k, 
                                var[*num_2d+v], *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    } else {
        redisReply *reply;
        for(j = *latids; j <= *latide; j+=*latstep) {
            sb_reset(key);
            sb_appendf(key, "hmset %s ", hkey);
            for(i = *lonids; i <= *lonide; i+=*lonstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.12f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.12f ", i, j, k, var[*num_2d+v], 
                                *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
            errno = redisClusterAppendCommand(cc, sb_concat(key));
            if (errno) {
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                errno = -1;
                goto err;
            }
            
        }

        for(j = *latids; j <= *latide; j+=*latstep) {
            int r = redisClusterGetReply(cc, (void **)&reply);
            if(r == REDIS_ERR) {
                printf("Redis Reply Error!\n");
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                errno = -1;
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            if (errno) {
                goto err;     
            }
        }
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_da_outputf(redisClusterContext *cc, char *hkey, char *varlist, 
        int *lonids, int *lonide, int *lonstep, 
        int *latids, int *latide, int *latstep, int *mpas_num_lev_start, 
        int *mpas_num_lev, int *zstep, int *num_2d, int *num_3d, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
	StringBuilder	*key = sb_create();
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j*n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    if (n * sec_i * sec_j < 518400) {
        sb_appendf(key, "hmset %s ", hkey);
        for(i = *lonids; i <= *lonide; i+=*lonstep) {
            for(j = *latids; j <= *latide; j+=*latstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.6f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.6f ", i, j, k, var[*num_2d+v], 
                                *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    } else {
        redisReply *reply;
        for(j = *latids; j <= *latide; j+=*latstep) {
            sb_reset(key);
            sb_appendf(key, "hmset %s ", hkey);
            for(i = *lonids; i <= *lonide; i+=*lonstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.6f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.6f ", i, j, k, var[*num_2d+v], 
                                *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
            errno = redisClusterAppendCommand(cc, sb_concat(key));
            if (errno) {
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                errno = -1;
                goto err;
            }

        }

        for(j = *latids; j <= *latide; j+=*latstep) {
            int r = redisClusterGetReply(cc, (void **)&reply);
            if(r == REDIS_ERR) {
                printf("Redis Reply Error!\n");
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                errno = -1;
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            if (errno) {
                goto err;     
            }
        }
    }
#endif
err:
    sb_free(key);
    return errno;
}



//-------------------------------------------------
// non-block interfaces
//-------------------------------------------------

int redis_hmsetf2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, int *its, 
        int *ite, int *istride, int *jts, int *jte, int *jstride, int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn + (j - *jts) * n + v)));
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }
    
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn + (j - *jts) * n + v)));
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:%d:0:%s %.6f ", i, j, var[v], 
                          *(buffer+((i - *its) * sec_jn + (j - *jts) * n + v)));
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmsetf_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s %.6f ", i, j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
    if (errno) {
        goto err;
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetf2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *nn, float *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            exit(-1);
        }
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:%d:0:%s ", i, j, var[v]); 
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnf(cc, reply, buffer, &mm);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetf_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            exit(-1);
        }
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnf(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s ", i, j, k, var[v]); 
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnf(cc, reply, buffer, &mm);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmsetd2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn
                          +(j - *jts) * n + v)));
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }
    
    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                          *(buffer+((i - *its) * sec_jn
                          +(j - *jts) * n + v)));
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
              for(v=0; v<n; v++) {
                  sb_appendf(key, "%d:%d:0:%s %.12f ", i, j, var[v], 
                          *(buffer+((i - *its) * sec_jn
                          +(j - *jts) * n + v)));
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmsetd_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmset %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s %.12f ", i, j, k, var[v], 
                            *(buffer+((i - *its) * sec_kjn
                            +(j - *jts) * sec_kn
                            +(k - *kms)*n + v)));
                }
            }
        }
    }

    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_check(cc, reply);
    freeReplyObject(reply);
    if (errno) {
        goto err;
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetd2d_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *nn, double *buffer)
{
    int i, j, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_jn = sec_j*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }

    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:0:%s ", j, var[v]); 
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(v=0; v<n; v++) {
                sb_appendf(key, "%d:%d:0:%s ", i, j, var[v]); 
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnd(cc, reply, buffer, &mm);
    freeReplyObject(reply);
    if (errno) {
        goto err;
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_hmgetd_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *its, int *ite, int *istride, int *jts, int *jte, int *jstride, 
        int *kms, int *kme, int *kstride, int *nn, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = *nn;
	StringBuilder	*key = sb_create();
    int sec_i = *ite - *its + 1;
    int sec_j = *jte - *jts + 1;
    int sec_k = *kme - *kms + 1;
    int sec_kjn = sec_k * sec_j*n;
    int sec_kn = sec_k*n;
    int size = (*ite - *its + 1) * (*jte - *jts + 1) * (*kme - *kms + 1) * n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_R
    redisReply *reply;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    int mm = 0;
    for(i = *its; i <= *ite; i+=*istride) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmget %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmget %s:%d ", hkey, i);
#endif
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%s ", j, k, var[v]); 
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_returnd(cc, reply, buffer, &mm);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    sb_appendf(key, "hmget %s ", hkey);
    for(i = *its; i <= *ite; i+=*istride) {
        for(j = *jts; j <= *jte; j+=*jstride) {
            for(k = *kms; k <= *kme; k+=*kstride) {
                for(v=0; v<n; v++) {
                    sb_appendf(key, "%d:%d:%d:%s ", i, j, k, var[v]); 
                }
            }
        }
    }

    int mm = 0;
    redisReply *reply = redisClusterCommand(cc, sb_concat(key));
    errno = reply_returnd(cc, reply, buffer, &mm);
    freeReplyObject(reply);
#endif
err:
    sb_free(key);
    return errno;
}

int redis_da_outputd_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *lonids, int *lonide, int *lonstep, 
        int *latids, int *latide, int *latstep, int *mpas_num_lev_start, 
        int *mpas_num_lev, int *zstep, int *num_2d, int *num_3d, double *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
	StringBuilder	*key = sb_create();
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j*n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }
#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.12f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.12f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    if (n * sec_i * sec_j < 518400) {
        sb_appendf(key, "hmset %s ", hkey);
        for(i = *lonids; i <= *lonide; i+=*lonstep) {
            for(j = *latids; j <= *latide; j+=*latstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.12f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.12f ", i, j, k, 
                                var[*num_2d+v], *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    } else {
        redisReply *reply;
        for(j = *latids; j <= *latide; j+=*latstep) {
            sb_reset(key);
            sb_appendf(key, "hmset %s ", hkey);
            for(i = *lonids; i <= *lonide; i+=*lonstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.12f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.12f ", i, j, k, var[*num_2d+v], 
                                *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
            errno = redisClusterAppendCommand(cc, sb_concat(key));
            if (errno) {
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                errno = -1;
                goto err;
            }
        }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

        for(j = *latids; j <= *latide; j+=*latstep) {
            int r = redisClusterGetReply(cc, (void **)&reply);
            if(r == REDIS_ERR) {
                printf("Redis Reply Error!\n");
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                errno = -1;
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            if (errno) {
                goto err;     
            }
        }
    }
#endif
err:
    sb_free(key);
    return errno;
}

int redis_da_outputf_nonblock(redisClusterContext *cc, char *hkey, char *varlist, 
        int *lonids, int *lonide, int *lonstep, 
        int *latids, int *latide, int *latstep, int *mpas_num_lev_start, 
        int *mpas_num_lev, int *zstep, int *num_2d, int *num_3d, float *buffer)
{
    int i, j, k, v;
    int errno = 0;
    int n = (*num_2d) + (*num_3d) * (*mpas_num_lev) ;
	StringBuilder	*key = sb_create();
    int sec_i = *lonide - *lonids + 1;
    int sec_j = *latide - *latids + 1;
    int sec_kjn = sec_j*n;

    /*string split*/
    char var[n][32];
    char* token = strtok( varlist, ":");
    i = 0;
    while( token != NULL )
    {
        strcpy(var[i], token); 
        token = strtok( NULL, ":");
        i++;
    }

#ifdef LON_SLICE
#ifdef REDIS_PIPE_W
    redisReply *reply;
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        errno = redisClusterAppendCommand(cc, sb_concat(key));
        if (errno) {
            err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
            errno = -1;
            goto err;
        }
    }

    errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        int r = redisClusterGetReply(cc, (void **)&reply);
        if(r == REDIS_ERR) {
            printf("Redis Reply Error!\n");
            err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
            errno = -1;
            freeReplyObject(reply);
            goto err;
        }
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;     
        }
    }
#else
    for(i = *lonids; i <= *lonide; i+=*lonstep) {
        sb_reset(key);
#ifdef LON_MAP
        sb_appendf(key, "hmset %s{%d} ", hkey, i);
#else
        sb_appendf(key, "hmset %s:%d ", hkey, i);
#endif
        for(j = *latids; j <= *latide; j+=*latstep) {
            for(v=0; v<*num_2d; v++) {
                sb_appendf(key, "%d:0:%s %.6f ", j, var[v], 
                        *(buffer+((i - *lonids) * sec_kjn
                        +(j - *latids) * n + v)));
            }
            for(v=0; v<*num_3d; v++) {
                for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                    sb_appendf(key, "%d:%d:%s %.6f ", j, k, var[*num_2d+v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n
                            +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                }
            }
        }
        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    }
#endif
#else
    if (n * sec_i * sec_j < 518400) {
        sb_appendf(key, "hmset %s ", hkey);
        for(i = *lonids; i <= *lonide; i+=*lonstep) {
            for(j = *latids; j <= *latide; j+=*latstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.6f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.6f ", i, j, k, var[*num_2d+v], 
                                *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
        }

        redisReply *reply = redisClusterCommand(cc, sb_concat(key));
        errno = reply_check(cc, reply);
        freeReplyObject(reply);
        if (errno) {
            goto err;
        }
    } else {
        redisReply *reply;
        for(j = *latids; j <= *latide; j+=*latstep) {
            sb_reset(key);
            sb_appendf(key, "hmset %s ", hkey);
            for(i = *lonids; i <= *lonide; i+=*lonstep) {
                for(v=0; v<*num_2d; v++) {
                    sb_appendf(key, "%d:%d:0:%s %.6f ", i, j, var[v], 
                            *(buffer+((i - *lonids) * sec_kjn
                            +(j - *latids) * n + v)));
                }
                for(v=0; v<*num_3d; v++) {
                    for(k = *mpas_num_lev_start; k <= *mpas_num_lev; k+=*zstep) {
                        sb_appendf(key, "%d:%d:%d:%s %.6f ", i, j, k, var[*num_2d+v], 
                                *(buffer+((i - *lonids) * sec_kjn
                                +(j - *latids) * n
                                +(*num_2d + v*(*mpas_num_lev)) + k - 1)));
                    }
                }
            }
            errno = redisClusterAppendCommand(cc, sb_concat(key));
            if (errno) {
                err_print("EXECUTE_COMMAND_ERROR:AppendCommandError:%s", cc->errstr);
                errno = -1;
                goto err;
            }
        }

	errno = redis_sendall(cc);
    if (errno) {
        err_print("EXECUTE_COMMAND_ERROR:SendAllError:%s", cc->errstr);
        goto err;
    }

        for(j = *latids; j <= *latide; j+=*latstep) {
            int r = redisClusterGetReply(cc, (void **)&reply);
            if(r == REDIS_ERR) {
                printf("Redis Reply Error!\n");
                err_print("EXECUTE_COMMAND_ERROR:GetReplyError:%s", cc->errstr);
                errno = -1;
                freeReplyObject(reply);
                goto err;
            }
            errno = reply_check(cc, reply);
            freeReplyObject(reply);
            if (errno) {
                goto err;     
            }
        }
    }
#endif
err:
    sb_free(key);
    return errno;
}

//-----------------------------------------------------------
// reply
//-----------------------------------------------------------
int reply_returnd(redisClusterContext *cc, redisReply *reply, double *buf, int *tag)
{
    int i;
	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}    switch(reply->type) {
        case 1 :
            *buf = atof(reply->str);
            break;
        case 2:
            for(i = 0; i < reply->elements; i++) {
                *(buf + *tag) = atof(reply->element[i]->str);
                (*tag)++;
            }
            break;
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;
    }
	return errno;
}

int reply_returnf(redisClusterContext *cc, redisReply *reply, float *buf, int *tag)
{
    int i;
 	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}   
    switch(reply->type) {
        case 1 :
            *buf = atof(reply->str);
            break;
        case 2:
            for(i = 0; i < reply->elements; i++) {
                *(buf + *tag) = atof(reply->element[i]->str);
                (*tag)++;
            }
            break;
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;    }
	return errno;
}

int reply_returni(redisClusterContext *cc, redisReply *reply, int *buf, int *tag)
{
    int i;
	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}
    switch(reply->type) {
        case 1 :
            *buf = atoi(reply->str);
            break;
        case 2:
            for(i = 0; i < reply->elements; i++) {
                *(buf + *tag) = atoi(reply->element[i]->str);
                (*tag)++;
            }
            break;
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;

    }
    return errno;
}

int reply_returns(redisClusterContext *cc, redisReply *reply, char *buf, int *len) {
	int i;
	int errno = 0;
	if (reply == NULL) {
        err_print("EXECUTE_COMMAND_ERROR:%s", cc->errstr);
		return -1;
	}
    switch(reply->type) {
        case 1 :
			memcpy(buf, reply->str, reply->len);
			*len = reply->len;
            break;
        case 2:
        case 3:
            printf("return integer\n");
            break;
        case 4:
            err_print("ERROR_FROM_REDIS: data is null! please check the database!");
			errno = -2;
            break;
        case 5:
            if(!(strcmp(reply->str,"OK")==0)) {
                err_print("ERROR_FROM_REDIS: return status: %s", reply->str);
				errno = -2;
            }
            break;
        case 6:
            err_print("ERROR_FROM_REDIS: %s!", reply->str);
			errno = -2;
            break;
        default:
            err_print("ERROR_FROM_REDIS: no match error please check redis data!\n");
			errno = -2;
    }
	return errno;
}

int reply_check(redisClusterContext *cc, redisReply *reply)
{
    if(NULL == reply) {
        err_print("Get Null Reply:%s", cc->errstr);
        freeReplyObject(reply);
		return -1;
    }
	return 0;
}
#endif
