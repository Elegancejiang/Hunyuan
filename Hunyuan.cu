
#include"include/HunyuanGRAPH.h"

/*Time function params*/
double part_all = 0;
struct timeval begin_part_all;
struct timeval   end_part_all;

double part_coarsen = 0;
struct timeval begin_part_coarsen;
struct timeval   end_part_coarsen;

double part_init = 0;
struct timeval begin_part_init;
struct timeval   end_part_init;

double part_uncoarsen = 0;
struct timeval begin_part_uncoarsen;
struct timeval   end_part_uncoarsen;

//four calculation pattern
double part_match = 0;
struct timeval begin_part_match;
struct timeval   end_part_match;

double part_contract = 0;
struct timeval begin_part_contract;
struct timeval   end_part_contract;
  


double part_cmatch = 0;
struct timeval begin_part_cmatch;
struct timeval   end_part_cmatch;

double part_ccontract = 0;
struct timeval begin_part_ccontract;
struct timeval   end_part_ccontract;

double part_bfs = 0;
struct timeval begin_part_bfs;
struct timeval   end_part_bfs;

double part_2refine = 0;
struct timeval begin_part_2refine;
struct timeval   end_part_2refine;

double part_2map = 0;
struct timeval begin_part_2map;
struct timeval   end_part_2map;

double part_slipt = 0;
struct timeval begin_part_slipt;
struct timeval   end_part_slipt;

double part_krefine = 0;
struct timeval begin_part_krefine;
struct timeval   end_part_krefine;

double part_map = 0;
struct timeval begin_part_map;
struct timeval   end_part_map;

double part_mallocrefine = 0;
struct timeval begin_part_mallocrefine;
struct timeval   end_part_mallocrefine;

//test
double krefine_atomicadd = 0;
struct timeval begin_krefine_atomicadd;
struct timeval   end_krefine_atomicadd;

//uncoarsen
double uncoarsen_Sum_maxmin_pwgts = 0;
double uncoarsen_Exnode_part1 = 0;
double uncoarsen_Exnode_part2 = 0;
struct timeval begin_general;
struct timeval   end_general;

double bndinfo_Find_real_bnd_info = 0;
double bndinfo_init_bnd_info = 0;
double bndinfo_find_kayparams = 0;
double bndinfo_initcucsr = 0;
double bndinfo_bb_segsort = 0;
double bndinfo_init_cu_que = 0;
double bndinfo_findcsr = 0;
struct timeval begin_bndinfo;
struct timeval   end_bndinfo;

//match
double sinitcuda_match = 0;
struct timeval begin_initcuda_match;
struct timeval   end_initcuda_match;

double scuda_match = 0;
struct timeval begin_cuda_match;
struct timeval   end_cuda_match;

double sfindc2 = 0;
struct timeval begin_findc2;
struct timeval   end_findc2;

double sinclusive_scan = 0;
struct timeval begin_inclusive_scan;
struct timeval   end_inclusive_scan;

double sfindc4 = 0;
struct timeval begin_findc4;
struct timeval   end_findc4;

//contract
double sexclusive_scan = 0;
struct timeval begin_exclusive_scan;
struct timeval   end_exclusive_scan;

double sfind_cnvtxsedge_original = 0;
struct timeval begin_find_cnvtxsedge_original;
struct timeval   end_find_cnvtxsedge_original;

double sbb_segsort = 0;
struct timeval begin_bb_segsort;
struct timeval   end_bb_segsort;

double sSort_cnedges_part1 = 0;
struct timeval begin_Sort_cnedges_part1;
struct timeval   end_Sort_cnedges_part1;

double sinclusive_scan2 = 0;
struct timeval begin_inclusive_scan2;
struct timeval   end_inclusive_scan2;

double sSort_cnedges_part2 = 0;
struct timeval begin_Sort_cnedges_part2;
struct timeval   end_Sort_cnedges_part2;

double sSort_cnedges_part2_5 = 0;
struct timeval begin_Sort_cnedges_part2_5;
struct timeval   end_Sort_cnedges_part2_5;

double sSort_cnedges_part3 = 0;
struct timeval begin_Sort_cnedges_part3;
struct timeval   end_Sort_cnedges_part3;

double coarsen_malloc = 0;
struct timeval begin_coarsen_malloc;
struct timeval   end_coarsen_malloc;

double coarsen_memcpy = 0;
struct timeval begin_coarsen_memcpy;
struct timeval   end_coarsen_memcpy;

double coarsen_free = 0;
struct timeval begin_coarsen_free;
struct timeval   end_coarsen_free;

double coarsen_else = 0;

double set_cpu_graph = 0;
struct timeval begin_set_cpu_graph;
struct timeval   end_set_cpu_graph;

double set_graph_1 = 0;
double set_graph_2 = 0;
double set_graph_3 = 0;
double set_graph_4 = 0;
struct timeval begin_set_graph;
struct timeval   end_set_graph;

double gpu_2way = 0;
struct timeval begin_gpu_2way;
struct timeval end_gpu_2way;

double malloc_2way = 0;
struct timeval begin_malloc_2way;
struct timeval end_malloc_2way;

double initmoveto = 0;
struct timeval begin_initmoveto;
struct timeval end_initmoveto;

double updatemoveto = 0;
struct timeval begin_updatemoveto;
struct timeval end_updatemoveto;

double computepwgts = 0;
struct timeval begin_computepwgts;
struct timeval end_computepwgts;

double thrustreduce = 0;
struct timeval begin_thrustreduce;
struct timeval end_thrustreduce;

double computegain = 0;
struct timeval begin_computegain;
struct timeval end_computegain;

double thrustsort = 0;
struct timeval begin_thrustsort;
struct timeval end_thrustsort;

double computegainv = 0;
struct timeval begin_computegainv;
struct timeval end_computegainv;

double inclusive = 0;
struct timeval begin_inclusive;
struct timeval end_inclusive;

double re_balance = 0;
struct timeval begin_rebalance;
struct timeval end_rebalance;

double malloc_split = 0;
struct timeval begin_malloc_split;
struct timeval end_malloc_split;

double memcpy_split = 0;
struct timeval begin_memcpy_split;
struct timeval end_memcpy_split;

double free_split = 0;
struct timeval begin_free_split;
struct timeval end_free_split;

double save_init = 0;
struct timeval begin_save_init;
struct timeval end_save_init;

/*Graph data structure*/
typedef struct hunyuangraph_graph_t {
	/*graph cpu params*/
	int nvtxs;                            //Graph vertex
	int nedges;	                          //Graph edge
	int *xadj;                            //Graph vertex csr array (xadj[nvtxs+1])
	int *adjncy;                          //Graph adjacency list (adjncy[nedges])
	int *adjwgt;   		                    //Graph edge weight array (adjwgt[nedges])
	int *vwgt;			                      //Graph vertex weight array(vwgr[nvtxs])
	int *tvwgt;                           //The sum of graph vertex weight 
	float *tvwgt_reverse;                 //The reciprocal of tvwgt
	int *label;                           //Graph vertex label(label[nvtxs])
	int *cmap;                            //The Label of graph vertex in cgraph(cmap[nvtxs]) 
	int mincut;                           //The min edfe-cut of graph partition
	int *where;                           //The label of graph vertex in which part(where[nvtxs]) 
	int *pwgts;                           //The partition vertex weight(pwgts[nparts])
	int nbnd;                             //Boundary vertex number
	int *bndlist;                         //Boundary vertex list
	int *bndptr;                          //Boundary vertex pointer
	int *id;                              //The sum of edge weight in same part
	int *ed;                              //The sum of edge weight in different part
	struct hunyuangraph_graph_t *coarser; //The coarser graph
	struct hunyuangraph_graph_t *finer;   //The finer graph
	/*graph gpu params*/
	int *cuda_xadj;
	int *cuda_adjncy;
	int *cuda_adjwgt;
	int *cuda_vwgt;               
	int *cuda_match;                      //CUDA graph vertex match array(match[nvtxs])
	int *cuda_cmap;
	// int *cuda_maxvwgt;                    //CUDA graph constraint vertex weight 
	int *txadj;                  //CUDA graph vertex pairs csr edge array(txadj[cnvtxs+1])
	//   int *cuda_real_nvtxs;                 //CUDA graph params (i<match[cmap[i]])
	//   int *cuda_s;                          //CUDA support array (cuda_s[nvtxs])
	int *tadjwgt;       //CUDA support scan array (tadjwgt[nedges])
	//   int *cuda_scan_nedges_original;       //CUDA support scan array (cuda_scan_nedges_original[nedges])
	int *tadjncy;      //CUDA support scan array (tadjncy[nedges])
	int *cuda_maxwgt;                     //CUDA part weight array (cuda_maxwgt[npart])
	int *cuda_minwgt;                     //CUDA part weight array (cuda_minwgt[npart])
	int *cuda_where;
	int *cuda_label;
	int *cuda_pwgts;
	int *cuda_bnd;
	int *cuda_bndnum;
	int *cpu_bndnum;
	int *cuda_info;                       //CUDA support array(cuda_info[bnd_num*nparts])
	int *cuda_real_bnd_num;
	int *cuda_real_bnd;
	//   int *cuda_tvwgt;  // graph->tvwgt
	float *cuda_tpwgts;

	/*Refinement available generate array*/
	int *cuda_bn;                             
	int *cuda_bt;
	// int *cuda_g;
	int *cuda_csr;
	int *cuda_que;
} hunyuangraph_graph_t;

/*Memory allocation information*/
typedef struct hunyuangraph_mop_t {
  int type;
  ssize_t nbytes;
  void *ptr;
} hunyuangraph_mop_t;

/*Algorithm information*/
typedef struct hunyuangraph_mcore_t {
  void *core;	
  size_t coresize;     
  size_t corecpos;            
  size_t nmops;         
  size_t cmop;         
  hunyuangraph_mop_t *mops;      
  size_t num_callocs;   
  size_t num_hallocs;   
  size_t size_callocs;  
  size_t size_hallocs;  
  size_t cur_callocs;   
  size_t cur_hallocs;  
  size_t max_callocs;   
  size_t max_hallocs;   

} hunyuangraph_mcore_t;

/*Control information*/
typedef struct hunyuangraph_admin_t {
  int Coarsen_threshold;		
  int nIparts;      
  int no2hop;                                                                                                                                 
  int iteration_num;                               
  int maxvwgt;		                
  int nparts;
  int ncuts;
  float *ubfactors;            
  float *tpwgts;               
  float *part_balance;               
  float cfactor;               
  hunyuangraph_mcore_t *mcore;    
  size_t nbrpoolsize;      
  size_t nbrpoolcpos;                  

} hunyuangraph_admin_t;

/*Heap information*/
typedef struct hunyuangraph_rkv_t{
  float key;
  int val;
} hunyuangraph_rkv_t;

/*Queue information*/
typedef struct {
  ssize_t nnodes;
  ssize_t maxnodes;
  hunyuangraph_rkv_t   *heap;
  ssize_t *locator;
} hunyuangraph_queue_t;

typedef struct {
  int key;
  int ptr;
} kp_t;

typedef struct ikv_t{
  int key;
  int val;
} ikv_t;

/*Define functions*/
#define ALIGNMENT 4
#define IDX_MAX   INT64_MAX
#define hunyuangraph_max(m,n) ((m)>=(n)?(m):(n))
#define hunyuangraph_min(m,n) ((m)>=(n)?(n):(m))
#define hunyuangraph_swap(m,n,temp) do{(temp)=(m);(m)=(n);(n)=(temp);} while(0) 
#define hunyuangraph_tocsr(i,n,c) do{for(i=1;i<n;i++)c[i]+= c[i-1];for(i=n;i>0;i--)c[i]=c[i-1];c[0]=0;} while(0)
#define SHIFTCSR(i, n, a) do {for (i=n; i>0; i--) a[i] = a[i-1]; a[0] = 0; } while(0) 
#define hunyuangraph_add_sub(m,n,temp) do{(m)+=(temp);(n)-=(temp);} while(0)
#define hunyuangraph_listinsert(n,list,lptr,i) do{list[n]=i;lptr[i]=(n)++;} while(0) 
#define hunyuangraph_listdelete(n,list,lptr,i) do{list[lptr[i]]=list[--(n)];lptr[list[n]]=lptr[i];lptr[i]=-1;} while(0) 
#define M_GT_N(m,n) ((m)>(n))

#define _GKQSORT_SWAP(a, b, t) ((void)((t = *a), (*a = *b), (*b = t)))
#define _GKQSORT_STACK_SIZE	    (8 * sizeof(size_t))
#define _GKQSORT_PUSH(top, low, high) (((top->_lo = (low)), (top->_hi = (high)), ++top))
#define	_GKQSORT_POP(low, high, top)  ((--top, (low = top->_lo), (high = top->_hi)))
#define	_GKQSORT_STACK_NOT_EMPTY	    (_stack < _top)

#define GK_MKQSORT(GKQSORT_TYPE,GKQSORT_BASE,GKQSORT_NELT,GKQSORT_LT)   \
{									\
  GKQSORT_TYPE *const _base = (GKQSORT_BASE);				\
  const size_t _elems = (GKQSORT_NELT);					\
  GKQSORT_TYPE _hold;							\
									\
  if (_elems == 0)                                                      \
    return;                                                             \
  if (_elems > 4) {					\
    GKQSORT_TYPE *_lo = _base;						\
    GKQSORT_TYPE *_hi = _lo + _elems - 1;				\
    struct {								\
      GKQSORT_TYPE *_hi; GKQSORT_TYPE *_lo;				\
    } _stack[_GKQSORT_STACK_SIZE], *_top = _stack + 1;			\
    while (_GKQSORT_STACK_NOT_EMPTY) {					\
      GKQSORT_TYPE *_left_ptr; GKQSORT_TYPE *_right_ptr;		\
      GKQSORT_TYPE *_mid = _lo + ((_hi - _lo) >> 1);			\
      if (GKQSORT_LT (_mid, _lo))					\
        _GKQSORT_SWAP (_mid, _lo, _hold);				\
      if (GKQSORT_LT (_hi, _mid))					\
        _GKQSORT_SWAP (_mid, _hi, _hold);				\
      else								\
        goto _jump_over;						\
      if (GKQSORT_LT (_mid, _lo))					\
        _GKQSORT_SWAP (_mid, _lo, _hold);				\
  _jump_over:;								\
      _left_ptr  = _lo + 1;						\
      _right_ptr = _hi - 1;						\
      do {								\
        while (GKQSORT_LT (_left_ptr, _mid))				\
         ++_left_ptr;							\
        while (GKQSORT_LT (_mid, _right_ptr))				\
          --_right_ptr;							\
        if (_left_ptr < _right_ptr) {					\
          _GKQSORT_SWAP (_left_ptr, _right_ptr, _hold);			\
          if (_mid == _left_ptr)					\
            _mid = _right_ptr;						\
          else if (_mid == _right_ptr)					\
            _mid = _left_ptr;						\
          ++_left_ptr;							\
          --_right_ptr;							\
        }								\
        else if (_left_ptr == _right_ptr) {				\
          ++_left_ptr;							\
          --_right_ptr;							\
          break;							\
        }								\
      } while (_left_ptr <= _right_ptr);				\
      if (_right_ptr - _lo <= 4) {			\
        if (_hi - _left_ptr <= 4)			\
          _GKQSORT_POP (_lo, _hi, _top);				\
        else								\
          _lo = _left_ptr;						\
      }									\
      else if (_hi - _left_ptr <= 4)			\
        _hi = _right_ptr;						\
      else if (_right_ptr - _lo > _hi - _left_ptr) {			\
        _GKQSORT_PUSH (_top, _lo, _right_ptr);				\
        _lo = _left_ptr;						\
      }									\
      else {								\
        _GKQSORT_PUSH (_top, _left_ptr, _hi);				\
        _hi = _right_ptr;						\
      }									\
    }									\
  }									\
  {									\
    GKQSORT_TYPE *const _end_ptr = _base + _elems - 1;			\
    GKQSORT_TYPE *_tmp_ptr = _base;					\
    register GKQSORT_TYPE *_run_ptr;					\
    GKQSORT_TYPE *_thresh;						\
    _thresh = _base + 4;				\
    if (_thresh > _end_ptr)						\
      _thresh = _end_ptr;						\
    for (_run_ptr = _tmp_ptr + 1; _run_ptr <= _thresh; ++_run_ptr)	\
      if (GKQSORT_LT (_run_ptr, _tmp_ptr))				\
        _tmp_ptr = _run_ptr;						\
    if (_tmp_ptr != _base)						\
      _GKQSORT_SWAP (_tmp_ptr, _base, _hold);				\
    _run_ptr = _base + 1;						\
    while (++_run_ptr <= _end_ptr) {					\
      _tmp_ptr = _run_ptr - 1;						\
      while (GKQSORT_LT (_run_ptr, _tmp_ptr))				\
        --_tmp_ptr;							\
      ++_tmp_ptr;							\
      if (_tmp_ptr != _run_ptr) {					\
        GKQSORT_TYPE *_trav = _run_ptr + 1;				\
        while (--_trav >= _run_ptr) {					\
          GKQSORT_TYPE *_hi; GKQSORT_TYPE *_lo;				\
          _hold = *_trav;						\
          for (_hi = _lo = _trav; --_lo >= _tmp_ptr; _hi = _lo)		\
            *_hi = *_lo;						\
          *_hi = _hold;							\
        }								\
      }									\
    }									\
  }									\
}

/*pointer*/
char *deviceMemory;
char *front_pointer;
char *back_pointer;
char *lmove_pointer;
char *rmove_pointer;
char *tmove_pointer;
size_t used_by_me_now = 0;
size_t used_by_me_max = 0;
size_t used_by_init_now = 0;
size_t used_by_init_max = 0;

void CPU_malloc(size_t size)
{
	used_by_init_now += size;
	used_by_init_max = hunyuangraph_max(used_by_init_now,used_by_init_max);
}

void CPU_free(size_t size)
{
	used_by_init_now -= size;
	used_by_init_max = hunyuangraph_max(used_by_init_now,used_by_init_max);
}

void Init_GPU_Memory(size_t remainingMem)
{
	front_pointer = deviceMemory;
	back_pointer  = deviceMemory + remainingMem;
	lmove_pointer = front_pointer;
	rmove_pointer = back_pointer;
}

void Malloc_GPU_Memory(size_t nvtxs, size_t nedges)
{
	// 获取设备属性
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);

	// 检测当前程序执行前显存已使用的大小
	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);
	size_t usableMem = totalMem - freeMem;

	// 计算剩余的部分显存大小
	// size_t remainingMem = freeMem - usableMem - 5 * nedges * sizeof(int);
	size_t remainingMem = freeMem - usableMem - (nvtxs + nedges) * sizeof(int);		//调试所用
	// size_t remainingMem = (freeMem - usableMem) / 2;	//调试所用
	// size_t remainingMem = 2 * 1024 * 1024 * 1024;	//调试所用
	if(remainingMem % 4 != 0) remainingMem += 4 - remainingMem % 4;

	// 分配对齐的显存空间(先不考虑对齐)
    // size_t pitch;
    // cudaMallocPitch(&deviceMemory, &pitch, remainingMem, 1);

	// 分配显存空间
	// size_t remainingMem = 2 * 1024 * 1024 * 1024;
	cudaMalloc(&deviceMemory, remainingMem);

	// 调整起始地址为对齐地址
    // size_t alignment = ALIGNMENT;  // 设置对齐要求，以字节为单位
    // front_pointer = (char*)deviceMemory;
    // front_pointer = (char*)(((size_t)front_pointer + alignment - 1) & ~(alignment - 1));

	Init_GPU_Memory(remainingMem);

	// printf("Total GPU Memory:                   %zuKB %zuMB %zuGB\n", totalMem / 1024,totalMem / 1024 / 1024,totalMem / 1024 / 1024 / 1024);
    // printf("Usable GPU Memory before execution: %zuKB %zuMB %zuGB\n", freeMem / 1024,freeMem / 1024 / 1024,freeMem / 1024 / 1024 / 1024);
    // printf("Malloc GPU Memory:                  %zuKB %zuMB %zuGB\n", remainingMem / 1024,remainingMem / 1024 / 1024,remainingMem / 1024 / 1024 / 1024);
    // printf("front_pointer=  %p\n", front_pointer);
    // printf("back_pointer=   %p\n", back_pointer);
}

void Free_GPU_Memory()
{
	// printf("lmove_pointer=  %p\n",lmove_pointer);
	// printf("rmove_pointer=  %p\n",rmove_pointer);
	// if(lmove_pointer == front_pointer) printf("The left stack has been freed\n");
	// else printf("error ------------ The left stack hasn't been freed\n");
	// if(rmove_pointer == back_pointer) printf("The right stack has been freed\n");
	// else printf("error ------------ The right stack hasn't been freed\n");
	printf("Max memory used of GPU: %10ldB %10ldKB %10ldMB %10ldGB\n",used_by_me_max,used_by_me_max / 1024,\
		used_by_me_max / 1024 / 1024, used_by_me_max/1024 / 1024 / 1024);
	cudaFree(deviceMemory);
}

// 申请和释放空间时参数为所需字节数
// l..为左端栈操作，r..为右端栈操作
void *lmalloc_with_check(size_t size, char *infor)
{
	if (size <= 0)
		printf("lmalloc_with_check( size <= 0)\n");

	void *malloc_address = NULL;

	tmove_pointer = lmove_pointer + size;
	// printf("lmove_pointer=%p\n",lmove_pointer);
	// printf("rmove_pointer=%p\n",rmove_pointer);
	// printf("size=         %d\n",size);
	// printf("tmove_pointer=%p\n",tmove_pointer);

	if (tmove_pointer > rmove_pointer)
		printf("error ------------ don't have enough GPU memory %s\n",infor);
	else
	{
		malloc_address = lmove_pointer;
		lmove_pointer = tmove_pointer;
		used_by_me_now += size;
		used_by_me_max = hunyuangraph_max(used_by_me_max,used_by_me_now);
	}

	// printf("lmalloc_with_check:%s\n",infor);
	// printf("used memory=    %10ld\n",used_by_me_now);
	// printf("malloc_address= %p\n",malloc_address);
	// printf("lmove_pointer=  %p\n",lmove_pointer);
	// printf("rmove_pointer=  %p lmalloc\n",rmove_pointer);
	// printf("available space %zuKB %zuMB %zuGB\n",(rmove_pointer - lmove_pointer) / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024 / 1024);
	// printf("\n");
	return malloc_address;
}

void *rmalloc_with_check(size_t size, char *infor)
{
	if (size <= 0)
		printf("rmalloc_with_check( size <= 0)\n");

	void *malloc_address = NULL;

	tmove_pointer = rmove_pointer - size;
	// printf("lmove_pointer=%p\n",lmove_pointer);
	// printf("rmove_pointer=%p\n",rmove_pointer);
	// printf("size=         %d\n",size);
	// printf("tmove_pointer=%p\n",tmove_pointer);

	if (tmove_pointer < lmove_pointer)
		printf("error ------------ don't have enough GPU memory %s\n",infor);
	else
	{
		malloc_address = tmove_pointer;
		rmove_pointer = tmove_pointer;
		used_by_me_now += size;
		used_by_me_max = hunyuangraph_max(used_by_me_max,used_by_me_now);
	}

	// printf("rmalloc_with_check:%s\n",infor);
	// printf("used memory=    %10ld\n",used_by_me_now);
	// printf("malloc_address= %p\n",malloc_address);
	// printf("lmove_pointer=  %p\n",lmove_pointer);
	// printf("rmove_pointer=  %p rmalloc\n",rmove_pointer);
	// printf("available space %zuKB %zuMB %zuGB\n",(rmove_pointer - lmove_pointer) / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024 / 1024);
	// printf("\n");
	return malloc_address;
}

void *lfree_with_check(size_t size, char *infor)
{
	if (size <= 0)
		printf("lfree_with_check( size <= 0)\n");

	void *malloc_address = NULL;

	tmove_pointer = lmove_pointer - size;

	// printf("lmove_pointer=%p\n",lmove_pointer);
	// printf("rmove_pointer=%p\n",rmove_pointer);
	// printf("size=         %d\n",size);
	// printf("tmove_pointer=%p\n",tmove_pointer);

	if (tmove_pointer < front_pointer)
		printf("error ------------ don't lmalloc enough GPU memory %s\n",infor);
	else
	{
		lmove_pointer = tmove_pointer;
		// malloc_address = lmove_pointer;
		used_by_me_now -= size;
		used_by_me_max = hunyuangraph_max(used_by_me_max,used_by_me_now);
	}

	// printf("lfree_with_check:%s\n",infor);
	// printf("used memory=    %10ld\n",used_by_me_now);
	// printf("lmove_pointer= %p\n",lmove_pointer);
	// printf("rmove_pointer= %p\n",rmove_pointer);
	// printf("malloc_address=lfree\n");
	// printf("available space %zuKB %zuMB %zuGB\n",(rmove_pointer - lmove_pointer) / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024 / 1024);
	// printf("\n");
	return malloc_address;
}

void *rfree_with_check(size_t size, char *infor)
{
	if (size <= 0)
		printf("rfree_with_check( size <= 0)\n");

	void *malloc_address = NULL;

	tmove_pointer = rmove_pointer + size;

	// printf("lmove_pointer=%p\n",lmove_pointer);
	// printf("rmove_pointer=%p\n",rmove_pointer);
	// printf("size=         %d\n",size);
	// printf("tmove_pointer=%p\n",tmove_pointer);

	if (tmove_pointer > back_pointer)
		printf("error ------------ don't rmalloc enough GPU memory %s\n",infor);
	else
	{
		rmove_pointer = tmove_pointer;
		// malloc_address = rmove_pointer;
		used_by_me_now -= size;
		used_by_me_max = hunyuangraph_max(used_by_me_max,used_by_me_now);
	}

	// printf("rfree_with_check:%s\n",infor);
	// printf("used memory=    %10ld\n",used_by_me_now);
	// printf("lmove_pointer=  %p\n",lmove_pointer);
	// printf("rmove_pointer=  %p\n",rmove_pointer);
	// printf("malloc_address= rfree\n");
	// printf("available space %zuKB %zuMB %zuGB\n",(rmove_pointer - lmove_pointer) / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024,(rmove_pointer - lmove_pointer) / 1024 / 1024 / 1024);
	// printf("\n");
	return malloc_address;
}

void ikvsorti(size_t n, ikv_t *base)
{
  #define ikey_lt(a, b) ((a)->key < (b)->key)
    GK_MKQSORT(ikv_t, base, n, ikey_lt);
  #undef ikey_lt
}

/*Compute log2 algorithm*/
int hunyuangraph_compute_log2(int a)
{
  int i;
  for(i=1;a>1;i++,a=a>>1);
  return i-1;
}

/*Get int rand number*/
int hunyuangraph_int_rand() 
{
  if(sizeof(int)<=sizeof(int32_t)) 
    return (int)(uint32_t)rand();
  else  
    return (int)(uint64_t)rand(); 
}

/*Get int rand number between (0,max)*/
int hunyuangraph_int_randinrange(int max) 
{
  return (int)((hunyuangraph_int_rand())%max); 
}

/*Compute sum of int array*/
int hunyuangraph_int_sum(size_t n, int *a)
{
  size_t i;
  int sum=0;
  for(i=0;i<n;i++,a+=1){
    sum+=(*a);
  }
  return sum;
}

/*Copy int array a to b*/
int  *hunyuangraph_int_copy(size_t n, int *a, int *b)
{
  return (int *)memmove((void *)b, (void *)a, sizeof(int)*n);
}

/*Set int array value*/
int *hunyuangraph_int_set_value(size_t n, int val, int *a)
{
  size_t i;
  for(i=0;i<n;i++){
    a[i]=val;
  }
  return a;
}

/*Compute sum of float array*/
float hunyuangraph_float_sum(size_t n, float *a)
{
  size_t i;
  float sum=0;
  for(i=0;i<n;i++,a+=1){
    sum+=(*a);
  }
  return sum;
}

/*Rescale tpwgts array*/
float *hunyuangraph_tpwgts_rescale(size_t n, float wsum, float *a)
{
  size_t i;
  for(i=0;i<n;i++,a+=1){
    (*a)*=wsum;
  }
  return a;
}

/*Compute Partition result edge-cut*/
int hunyuangraph_computecut(hunyuangraph_graph_t *graph, int *where)
{
  	int i,j,cut=0;
    for(i=0;i<graph->nvtxs;i++)
	{
		// printf("i=%d\n",i);
      	for(j=graph->xadj[i];j<graph->xadj[i+1];j++)
			if(where[i]!=where[graph->adjncy[j]])
				cut+=graph->adjwgt[j];
    }
  	return cut/2;
}

/*Set graph admin params*/
hunyuangraph_admin_t *hunyuangraph_set_graph_admin(int nparts, float *tpwgts, float *ubvec)
{
  int i;
  hunyuangraph_admin_t *hunyuangraph_admin;
  hunyuangraph_admin=(hunyuangraph_admin_t *)malloc(sizeof(hunyuangraph_admin_t));
  memset((void *)hunyuangraph_admin,0,sizeof(hunyuangraph_admin_t));

  hunyuangraph_admin->iteration_num=10;
  hunyuangraph_admin->Coarsen_threshold=200;
  hunyuangraph_admin->nparts=nparts; 

  hunyuangraph_admin->maxvwgt=0;  
  hunyuangraph_admin->ncuts=1; 

  hunyuangraph_admin->tpwgts=(float*)malloc(sizeof(float)*nparts);
  for(i=0;i<nparts;i++){
    hunyuangraph_admin->tpwgts[i]=1.0/nparts;
  }

  hunyuangraph_admin->ubfactors=(float*)malloc(sizeof(float));
  hunyuangraph_admin->ubfactors[0] =1.03;

  hunyuangraph_admin->part_balance =(float*) malloc(sizeof(float)*nparts);
  return hunyuangraph_admin;  
}

/*Set graph params*/
void hunyuangraph_init_cpu_graph(hunyuangraph_graph_t *graph) 
{
  memset((void *)graph,0,sizeof(hunyuangraph_graph_t));
  graph->nvtxs     = -1;
  graph->nedges    = -1;
  graph->xadj      = NULL;
  graph->vwgt      = NULL;
  graph->adjncy    = NULL;
  graph->adjwgt    = NULL;
  graph->label     = NULL;
  graph->cmap      = NULL;
  graph->tvwgt     = NULL;
  graph->tvwgt_reverse  = NULL;
  graph->where     = NULL;
  graph->pwgts     = NULL;
  graph->mincut    = -1;
  graph->nbnd      = -1;
  graph->id        = NULL;
  graph->ed        = NULL;
  graph->bndptr    = NULL;
  graph->bndlist   = NULL;
  graph->coarser   = NULL;
  graph->finer     = NULL;
}

/*Malloc graph*/
hunyuangraph_graph_t *hunyuangraph_create_cpu_graph(void)
{
    hunyuangraph_graph_t *graph = (hunyuangraph_graph_t *)malloc(sizeof(hunyuangraph_graph_t));
    hunyuangraph_init_cpu_graph(graph);
    return graph;
}

/*Set graph tvwgt value*/
void hunyuangraph_set_graph_tvwgt(hunyuangraph_graph_t *graph)
{
    if(graph->tvwgt == NULL)
	{
        graph->tvwgt = (int*)malloc(sizeof(int));
		CPU_malloc(sizeof(int));
	}

    if(graph->tvwgt_reverse == NULL)
	{
        graph->tvwgt_reverse = (float*)malloc(sizeof(float));
		CPU_malloc(sizeof(float));
	}

    graph->tvwgt[0] = hunyuangraph_int_sum(graph->nvtxs,graph->vwgt);
    graph->tvwgt_reverse[0] = 1.0 / (graph->tvwgt[0] > 0 ? graph->tvwgt[0] : 1);
}

/*Set graph vertex label*/
void hunyuangraph_set_graph_label(hunyuangraph_graph_t *graph)
{
    if(graph->label == NULL)
	{
        graph->label = (int *)malloc(sizeof(int) * (graph->nvtxs));
		CPU_malloc(sizeof(int) * graph->nvtxs);
	}

    for(int i = 0;i < graph->nvtxs;i++)
        graph->label[i] = i;
}

/*Set graph information*/
hunyuangraph_graph_t *hunyuangraph_set_graph(hunyuangraph_admin_t *hunyuangraph_admin, int nvtxs, int *xadj, int *adjncy, int *vwgt , int *adjwgt, int *tvwgt) 
{
    hunyuangraph_graph_t *graph = hunyuangraph_create_cpu_graph();

    graph->nvtxs  = nvtxs;
    graph->nedges = xadj[nvtxs];
    graph->xadj   = xadj;
    graph->adjncy = adjncy;

    graph->vwgt   = vwgt;
    graph->adjwgt = adjwgt;

    graph->tvwgt         = (int*)malloc(sizeof(int));
    graph->tvwgt_reverse = (float*)malloc(sizeof(float));

    graph->tvwgt[0]         = tvwgt[0];
    graph->tvwgt_reverse[0] = 1.0 / (graph->tvwgt[0] > 0 ? graph->tvwgt[0] : 1);

    //init label spend much time
    // graph->label = (int *)malloc(sizeof(int) * (graph->nvtxs));
    // for(int i = 0;i < graph->nvtxs;i++)
    //     graph->label[i] = i;
  
  	return graph;
}

hunyuangraph_graph_t *hunyuangraph_set_first_level_graph(int nvtxs, int *xadj, int *adjncy, int *vwgt , int *adjwgt)
{
    int i;
    hunyuangraph_graph_t *graph;
    
    graph = hunyuangraph_create_cpu_graph();
    graph->nvtxs  = nvtxs;
    graph->nedges = xadj[nvtxs];
    graph->xadj   = xadj;
    graph->adjncy = adjncy;
    
    graph->vwgt   = vwgt;
    graph->adjwgt = adjwgt;
    
    graph->tvwgt         = (int*)malloc(sizeof(int));
    graph->tvwgt_reverse = (float*)malloc(sizeof(float));
    graph->tvwgt[0]         = nvtxs;
    graph->tvwgt_reverse[0] = 1.0 / (graph->tvwgt[0] > 0 ? graph->tvwgt[0]:1);
    
    return graph;
} 

/*Creates mcore*/
hunyuangraph_mcore_t *hunyuangraph_create_mcore(size_t coresize)
{
  hunyuangraph_mcore_t *mcore;
  mcore=(hunyuangraph_mcore_t *)malloc(sizeof(hunyuangraph_mcore_t));
  memset(mcore,0,sizeof(hunyuangraph_mcore_t));

  mcore->coresize=coresize;
  mcore->corecpos=0;
  mcore->core=(coresize==0?NULL:(size_t*)malloc(sizeof(size_t)*(mcore->coresize)));
  mcore->nmops=2048;
  mcore->cmop=0;
  mcore->mops=(hunyuangraph_mop_t *)malloc((mcore->nmops)*sizeof(hunyuangraph_mop_t));

  return mcore;
}

/*Allocate work space*/
void hunyuangraph_allocatespace(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
  size_t coresize;
  coresize=3*(graph->nvtxs+1)*sizeof(int)+5*(hunyuangraph_admin->nparts+1)*sizeof(int)\
  +5*(hunyuangraph_admin->nparts+1)*sizeof(float);

  hunyuangraph_admin->mcore=hunyuangraph_create_mcore(coresize);
  hunyuangraph_admin->nbrpoolsize=0;
  hunyuangraph_admin->nbrpoolcpos=0;
}

/*Add memory allocation*/
void hunyuangraph_add_mcore(hunyuangraph_mcore_t *mcore, int type, size_t nbytes, void *ptr)
{
  if(mcore->cmop==mcore->nmops){
    mcore->nmops*=2;
    mcore->mops=(hunyuangraph_mop_t*)realloc(mcore->mops, mcore->nmops*sizeof(hunyuangraph_mop_t));
    if(mcore->mops==NULL){
      exit(0);
    }
  }

  mcore->mops[mcore->cmop].type=type;
  mcore->mops[mcore->cmop].nbytes=nbytes;
  mcore->mops[mcore->cmop].ptr=ptr;
  mcore->cmop++;

  switch(type){
    case 1:
      break;
    
    case 2:
      mcore->num_callocs++;
      mcore->size_callocs+=nbytes;
      mcore->cur_callocs+=nbytes;
      if(mcore->max_callocs<mcore->cur_callocs){
        mcore->max_callocs=mcore->cur_callocs;
      }
      break;
    
    case 3:
      mcore->num_hallocs++;
      mcore->size_hallocs+=nbytes;
      mcore->cur_hallocs+=nbytes;
      if(mcore->max_hallocs<mcore->cur_hallocs){
        mcore->max_hallocs=mcore->cur_hallocs;
      }
      break;
    
    default:
      exit(0);
  }
}

/*Malloc mcore*/
void *hunyuangraph_malloc_mcore(hunyuangraph_mcore_t *mcore, size_t nbytes)
{
  void *ptr;
  nbytes+=(nbytes%8==0?0:8-nbytes%8);

  if(mcore->corecpos+nbytes<mcore->coresize){
    ptr=((char *)mcore->core)+mcore->corecpos;
    mcore->corecpos+=nbytes;
    hunyuangraph_add_mcore(mcore,2,nbytes,ptr);
  }
  else{
    ptr=(size_t*)malloc(nbytes);
    hunyuangraph_add_mcore(mcore,3,nbytes,ptr);
  }

  return ptr;
}

/*Malloc mcore space*/
void *hunyuangraph_malloc_space(hunyuangraph_admin_t *hunyuangraph_admin, size_t nbytes)
{
  return hunyuangraph_malloc_mcore(hunyuangraph_admin->mcore,nbytes);
}

/*Malloc int mcore space*/
int *hunyuangraph_int_malloc_space(hunyuangraph_admin_t *hunyuangraph_admin, size_t n)
{
  return (int *)hunyuangraph_malloc_space(hunyuangraph_admin, n*sizeof(int));
}

ikv_t *ikvwspacemalloc(hunyuangraph_admin_t *hunyuangraph_admin, size_t nunmatched)
{
	return (ikv_t *)hunyuangraph_malloc_space(hunyuangraph_admin, nunmatched*sizeof(ikv_t));
}

/*Malloc float mcore space*/
float *hunyuangraph_float_malloc_space(hunyuangraph_admin_t *hunyuangraph_admin)
{
  return (float *)hunyuangraph_malloc_space(hunyuangraph_admin,2*sizeof(float));
}

/*Compute 2way balance params*/
void hunyuangraph_compute_2way_balance(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, float *tpwgts)
{
  int i;
  for(i=0;i<2;i++){
      hunyuangraph_admin->part_balance[i]=graph->tvwgt_reverse[0]/tpwgts[i];
  }
}

/*Get random permute of p*/
void hunyuangraph_int_randarrayofp(int n, int *p, int m, int flag)
{
  int i,u,v;
  int temp;
  if(flag==1){
    for(i=0;i<n;i++)
      p[i] = (int)i;
  }

  if(n<10){
    for(i=0;i<n;i++){

      v=hunyuangraph_int_randinrange(n);
      u=hunyuangraph_int_randinrange(n);
     
      hunyuangraph_swap(p[v],p[u],temp);

    }
  }
  else{
    for(i=0;i<m;i++){

      v=hunyuangraph_int_randinrange(n-3);
      u=hunyuangraph_int_randinrange(n-3);
      
      hunyuangraph_swap(p[v+0],p[u+2],temp);
      hunyuangraph_swap(p[v+1],p[u+3],temp);
      hunyuangraph_swap(p[v+2],p[u+0],temp);
      hunyuangraph_swap(p[v+3],p[u+1],temp);

    }
  }
}

/*Get permutation array*/
void hunyuangraph_matching_sort(hunyuangraph_admin_t *hunyuangraph_admin, int n, \
int max, int *keys, int *tperm, int *perm)
{
  int i,ii;
  int *counts;
  counts=hunyuangraph_int_set_value(max+2,0,hunyuangraph_int_malloc_space(hunyuangraph_admin,max+2));
  CPU_malloc(sizeof(int) * (max + 2));
  
  for(i=0; i<n; i++){
    counts[keys[i]]++;
  }
  
  hunyuangraph_tocsr(i,max+1,counts);
  
  for(ii=0;ii<n;ii++){
    i=tperm[ii];
    perm[counts[keys[i]]++]=i;
  }
}

/*Malloc cpu coarsen graph params*/
hunyuangraph_graph_t *hunyuangraph_set_cpu_cgraph(hunyuangraph_graph_t *graph, int cnvtxs)
{
  hunyuangraph_graph_t *cgraph;
  cgraph=hunyuangraph_create_cpu_graph();
  
  cgraph->nvtxs=cnvtxs;
  cgraph->xadj=(int*)malloc(sizeof(int)*(cnvtxs+1));
  cgraph->adjncy=(int*)malloc(sizeof(int)*(graph->nedges));
  cgraph->adjwgt=(int*)malloc(sizeof(int)*(graph->nedges));
  cgraph->vwgt=(int*)malloc(sizeof(int)*cnvtxs);
  cgraph->tvwgt=(int*)malloc(sizeof(int));
  cgraph->tvwgt_reverse=(float*)malloc(sizeof(float)); 
  
  cgraph->finer=graph;
  graph->coarser=cgraph;
  
  return cgraph;
}

/*Malloc gpu coarsen graph params*/
hunyuangraph_graph_t *hunyuangraph_set_gpu_cgraph(hunyuangraph_graph_t *graph, int cnvtxs)
{
	hunyuangraph_graph_t *cgraph = hunyuangraph_create_cpu_graph();
	
	cgraph->nvtxs=cnvtxs;
	//   cgraph->xadj=(int*)malloc(sizeof(int)*(cnvtxs+1));
	cgraph->tvwgt=(int*)malloc(sizeof(int));
	cgraph->tvwgt_reverse=(float*)malloc(sizeof(float)); 
	
	cgraph->finer=graph;
	graph->coarser=cgraph;
	
	return cgraph;
}

/*Create cpu coarsen graph by contract*/
void hunyuangraph_cpu_create_cgraph(hunyuangraph_admin_t *hunyuangraph_admin, \
hunyuangraph_graph_t *graph, int cnvtxs, int *match)
{
  int j,k,m,istart,iend,nvtxs,nedges,cnedges,v,u;
  int *xadj,*vwgt,*adjncy,*adjwgt;
  int *cmap,*htable;
  int *cxadj,*cvwgt,*cadjncy,*cadjwgt;
  hunyuangraph_graph_t *cgraph;
  
  nvtxs=graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
  cmap=graph->cmap;                  
  
  cgraph=hunyuangraph_set_cpu_cgraph(graph,cnvtxs);
  CPU_malloc(sizeof(int) * (cnvtxs + 1 + graph->nedges + graph->nedges + cnvtxs + 1));
  cxadj=cgraph->xadj;
  cvwgt=cgraph->vwgt;
  cadjncy=cgraph->adjncy;
  cadjwgt=cgraph->adjwgt;                               
  htable=hunyuangraph_int_set_value(cnvtxs,-1,hunyuangraph_int_malloc_space(hunyuangraph_admin,cnvtxs));
  CPU_malloc(sizeof(int) * cnvtxs);      
  cxadj[0] = cnvtxs = cnedges = 0; 
  nedges=graph->nedges;
   
  for(v=0;v<nvtxs;v++){

    if((u=match[v])<v)         
      continue;   

    cvwgt[cnvtxs]=vwgt[v];                 
    nedges=0;                                                    
    istart=xadj[v];
    iend=xadj[v+1];    

    for(j=istart;j<iend;j++){

      k=cmap[adjncy[j]];     

      if((m=htable[k])==-1){
        cadjncy[nedges]=k;                           
        cadjwgt[nedges] = adjwgt[j];                      
        htable[k] = nedges++;  
      }
      else{
        cadjwgt[m] += adjwgt[j];                                 
      }
    }

    if(v!=u){ 
      cvwgt[cnvtxs]+=vwgt[u];                   
      istart=xadj[u];                                    
      iend=xadj[u+1];      

      for(j=istart;j<iend;j++){
        k=cmap[adjncy[j]];

        if((m=htable[k])==-1){
          cadjncy[nedges]=k;
          cadjwgt[nedges]=adjwgt[j];
          htable[k]=nedges++;
        }
        else{
          cadjwgt[m] += adjwgt[j];
        }
      }

      if((j=htable[cnvtxs])!=-1){
        cadjncy[j]=cadjncy[--nedges];
        cadjwgt[j]=cadjwgt[nedges];
        htable[cnvtxs] = -1;
      }
    }

    for(j=0;j<nedges;j++){
       htable[cadjncy[j]] = -1;  
    }

    cnedges+=nedges;
    cxadj[++cnvtxs]=cnedges;
    cadjncy+=nedges;                                                                 
    cadjwgt+=nedges;
  }

  cgraph->nedges=cnedges;
  cgraph->tvwgt[0]=hunyuangraph_int_sum(cgraph->nvtxs,cgraph->vwgt); 
  cgraph->tvwgt_reverse[0]=1.0/(cgraph->tvwgt[0]>0?cgraph->tvwgt[0]:1);    

//   printf("cnvtxs=%d cnedges=%d\n",cgraph->nvtxs,cgraph->nedges);

}

/*************************************************************************/
/*! This function matches the unmatched vertices whose degree is less than
    maxdegree using a 2-hop matching that involves vertices that are two 
    hops away from each other. 
    The requirement of the 2-hop matching is a simple non-empty overlap
    between the adjancency lists of the vertices. */
/**************************************************************************/
int Match_2HopAny(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, int *perm, int *match, 
          int cnvtxs, size_t *r_nunmatched, size_t maxdegree)
{
  int i, pi, ii, j, jj, k, nvtxs;
  int *xadj, *adjncy, *colptr, *rowind;
  int *cmap;
  size_t nunmatched;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  cmap   = graph->cmap;

  nunmatched = *r_nunmatched;

  /* create the inverted index */
  colptr = hunyuangraph_int_set_value(nvtxs, 0, hunyuangraph_int_malloc_space(hunyuangraph_admin, nvtxs+1));
  CPU_malloc(sizeof(int) * (nvtxs + 1));
  for (i=0; i<nvtxs; i++) {
    if (match[i] == -1 && xadj[i+1]-xadj[i] < maxdegree) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        colptr[adjncy[j]]++;
    }
  }
  hunyuangraph_tocsr(i, nvtxs, colptr);

  rowind = hunyuangraph_int_malloc_space(hunyuangraph_admin, colptr[nvtxs]);
  CPU_malloc(sizeof(int) * (colptr[nvtxs]));
  for (pi=0; pi<nvtxs; pi++) {
    i = perm[pi];
    if (match[i] == -1 && xadj[i+1]-xadj[i] < maxdegree) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        rowind[colptr[adjncy[j]]++] = i;
    }
  }
  SHIFTCSR(i, nvtxs, colptr);

  /* compute matchings by going down the inverted index */
  for (pi=0; pi<nvtxs; pi++) {
    i = perm[pi];
    if (colptr[i+1]-colptr[i] < 2)
      continue;

    for (jj=colptr[i+1], j=colptr[i]; j<jj; j++) {
      if (match[rowind[j]] == -1) {
        for (jj--; jj>j; jj--) {
          if (match[rowind[jj]] == -1) {
            cmap[rowind[j]] = cmap[rowind[jj]] = cnvtxs++;
            match[rowind[j]]  = rowind[jj];
            match[rowind[jj]] = rowind[j];
            nunmatched -= 2;
            break;
          }
        }
      }
    }
  }

  *r_nunmatched = nunmatched;

  return cnvtxs;
}

/*************************************************************************/
/*! This function matches the unmatched vertices whose degree is less than
    maxdegree using a 2-hop matching that involves vertices that are two 
    hops away from each other. 
    The requirement of the 2-hop matching is that of identical adjacency
    lists.
 */
/**************************************************************************/
int Match_2HopAll(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, int *perm, int *match, 
          int cnvtxs, size_t *r_nunmatched, size_t maxdegree)
{
  int i, pi, pk, ii, j, jj, k, nvtxs, mask, idegree;
  int *xadj, *adjncy;
  int *cmap, *mark;
  ikv_t *keys;
  size_t nunmatched, ncand;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  cmap   = graph->cmap;

  nunmatched = *r_nunmatched;
  mask = IDX_MAX/maxdegree;

  /* collapse vertices with identical adjancency lists */
  keys = ikvwspacemalloc(hunyuangraph_admin, nunmatched);
  CPU_malloc(sizeof(ikv_t) * nunmatched);
  for (ncand=0, pi=0; pi<nvtxs; pi++) {
    i = perm[pi];
    idegree = xadj[i+1]-xadj[i];
    if (match[i] == -1 && idegree > 1 && idegree < maxdegree) {
      for (k=0, j=xadj[i]; j<xadj[i+1]; j++) 
        k += adjncy[j]%mask;
      keys[ncand].val = i;
      keys[ncand].key = (k%mask)*maxdegree + idegree;
      ncand++;
    }
  }
  ikvsorti(ncand, keys);

  match=hunyuangraph_int_set_value(nvtxs,0, hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs));
  CPU_malloc(sizeof(int) * nvtxs);
  for (pi=0; pi<ncand; pi++) {
    i = keys[pi].val;
    if (match[i] != -1)
      continue;

    for (j=xadj[i]; j<xadj[i+1]; j++)
      mark[adjncy[j]] = i;

    for (pk=pi+1; pk<ncand; pk++) {
      k = keys[pk].val;
      if (match[k] != -1)
        continue;

      if (keys[pi].key != keys[pk].key)
        break;
      if (xadj[i+1]-xadj[i] != xadj[k+1]-xadj[k])
        break;

      for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
        if (mark[adjncy[jj]] != i)
          break;
      }
      if (jj == xadj[k+1]) {
        cmap[i] = cmap[k] = cnvtxs++;
        match[i] = k;
        match[k] = i;
        nunmatched -= 2;
        break;
      }
    }
  }

  *r_nunmatched = nunmatched;

  return cnvtxs;
}

/*************************************************************************/
/*! This function matches the unmatched vertices using a 2-hop matching 
    that involves vertices that are two hops away from each other. */
/**************************************************************************/
int Match_2Hop(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, int *perm, int *match, 
          int cnvtxs, size_t nunmatched)
{

	cnvtxs = Match_2HopAny(hunyuangraph_admin, graph, perm, match, cnvtxs, &nunmatched, 2);
	cnvtxs = Match_2HopAll(hunyuangraph_admin, graph, perm, match, cnvtxs, &nunmatched, 64);
	if (nunmatched > 1.5*0.1*graph->nvtxs) 
		cnvtxs = Match_2HopAny(hunyuangraph_admin, graph, perm, match, cnvtxs, &nunmatched, 3);
	if (nunmatched > 2.0*0.1*graph->nvtxs) 
		cnvtxs = Match_2HopAny(hunyuangraph_admin, graph, perm, match, cnvtxs, &nunmatched, graph->nvtxs);

 	return cnvtxs;
}

int hunyuangraph_cpu_match_RM(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
	int i, pi, ii, j, jj, jjinc, k, nvtxs, cnvtxs, maxidx, last_unmatched;
	int *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
	int *match, *cmap, *perm;
	size_t nunmatched=0;

	nvtxs  = graph->nvtxs;
	xadj   = graph->xadj;
	vwgt   = graph->vwgt;
	adjncy = graph->adjncy;
	adjwgt = graph->adjwgt;
	cmap   = graph->cmap;

	maxvwgt = (int *)malloc(sizeof(int));
	maxvwgt[0] = hunyuangraph_admin->maxvwgt;

	match=hunyuangraph_int_set_value(nvtxs,-1, hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs));
  	perm=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
	CPU_malloc(sizeof(int) * nvtxs * 2);

	hunyuangraph_int_randarrayofp(nvtxs,perm,nvtxs/8,1);   

	for (cnvtxs=0, last_unmatched=0, pi=0; pi<nvtxs; pi++) 
	{
		i = perm[pi];

		if (match[i] == -1) {  /* Unmatched */
			maxidx = i;

			if (vwgt[i] < maxvwgt[0]) {
				/* Deal with island vertices. Find a non-island and match it with. 
				The matching ignores ctrl->maxvwgt requirements */
				if (xadj[i] == xadj[i+1]) {
					last_unmatched = hunyuangraph_max(pi, last_unmatched)+1;
					for (; last_unmatched<nvtxs; last_unmatched++) {
						j = perm[last_unmatched];
						if (match[j] == -1) {
							maxidx = j;
							break;
						}
					}
				}
				else {
				/* Find a random matching, subject to maxvwgt constraints */
					/* single constraint version */
					for (j=xadj[i]; j<xadj[i+1]; j++) {
						k = adjncy[j];
						if (match[k] == -1 && vwgt[i]+vwgt[k] <= maxvwgt[0]) {
							maxidx = k;
							break;
						}
					}

					/* If it did not match, record for a 2-hop matching. */
					if (maxidx == i && 3*vwgt[i] < maxvwgt[0]) {
						nunmatched++;
						maxidx = -1;
					}
				}
			}

			if (maxidx != -1) {
				cmap[i]  = cmap[maxidx] = cnvtxs++;
				match[i] = maxidx;
				match[maxidx] = i;
			}
		}
	}

	/* see if a 2-hop matching is required/allowed */
	if (!hunyuangraph_admin->no2hop && nunmatched > 0.1*nvtxs) {
		cnvtxs = Match_2Hop(hunyuangraph_admin, graph, perm, match, cnvtxs, nunmatched);
	}

	/* match the final unmatched vertices with themselves and reorder the vertices 
     of the coarse graph for memory-friendly contraction */
	for (cnvtxs=0, i=0; i<nvtxs; i++) {
		if (match[i] == -1) {
			match[i] = i;
			cmap[i]  = cnvtxs++;
		}
		else {
			if (i <= match[i]) 
				cmap[i] = cmap[match[i]] = cnvtxs++;
		}
	}

	hunyuangraph_cpu_create_cgraph(hunyuangraph_admin, graph, cnvtxs, match);
}

/*Get cpu graph matching params by hem*/
int hunyuangraph_cpu_match_HEM(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
  cudaDeviceSynchronize();
  gettimeofday(&begin_part_cmatch,NULL);

  int i,j,pi,k,nvtxs,cnvtxs,maxidx,maxwgt,last_unmatched,aved;
  int *xadj,*vwgt,*adjncy,*adjwgt,maxvwgt;
  int *match,*cmap,*d,*perm,*tperm;
  size_t nunmatched=0;

  nvtxs=graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
  cmap=graph->cmap;
  maxvwgt=hunyuangraph_admin->maxvwgt;
  
  cnvtxs=0;
  match=hunyuangraph_int_set_value(nvtxs,-1, hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs));
  perm=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  tperm=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  d=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);         
  CPU_malloc(sizeof(int) * nvtxs * 4);
  hunyuangraph_int_randarrayofp(nvtxs,tperm,nvtxs/8,1);   
  aved=0.7*(xadj[nvtxs]/nvtxs);

  for(i=0;i<nvtxs;i++){ 
    d[i]=(xadj[i+1]-xadj[i]>aved?aved:xadj[i+1]-xadj[i]);
  }

  hunyuangraph_matching_sort(hunyuangraph_admin,nvtxs,aved,d,tperm,perm);         
  
  last_unmatched=0;
  for(pi=0;pi<nvtxs;pi++) 
  {
    i=perm[pi];  

    if(match[i]==-1){  
      maxidx=i;                                                                               
      maxwgt=-1;           

	  if(vwgt[i] < maxvwgt)
	  {
		/* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) { 
			last_unmatched = hunyuangraph_max(pi, last_unmatched)+1;
			for (; last_unmatched<nvtxs; last_unmatched++) {
				j = perm[last_unmatched];
				if (match[j] == -1) {
					maxidx = j;
					break;
				}
			}
        }
		else
		{
			/* Find a heavy-edge matching, subject to maxvwgt constraints */
			/* single constraint version */
			for(j=xadj[i];j<xadj[i+1];j++){
				k=adjncy[j];

				if(match[k]==-1&&maxwgt<adjwgt[j]&&vwgt[i]+vwgt[k]<=maxvwgt){
					maxidx=k;
					maxwgt=adjwgt[j];
				}   
			}

			if(maxidx==i&&3*vwgt[i]<maxvwgt){ 
				maxidx=-1;
			}
	    }
	  }

      if(maxidx!=-1){
        cmap[i]=cmap[maxidx]=cnvtxs++;              
        match[i]=maxidx;                                        
        match[maxidx]=i; 
      }
    }
  }

  /* see if a 2-hop matching is required/allowed */
	if (!hunyuangraph_admin->no2hop && nunmatched > 0.1*nvtxs) 
	{
		cnvtxs = Match_2Hop(hunyuangraph_admin, graph, perm, match, cnvtxs, nunmatched);
	}

  for(cnvtxs=0,i=0;i<nvtxs;i++){
    if(match[i]==-1){
      match[i]=i;
      cmap[i]=cnvtxs++;                                                    
    }
    else{
      if(i<=match[i]){ 
        cmap[i]=cmap[match[i]]=cnvtxs++;
      }
    }
  }

  /*for(i = 0;i < nvtxs;i += 2)
  {
	if(i + 1 < nvtxs) 
	{
		match[i] = i + 1;
		match[i + 1] = i;
		cmap[i] = i / 2;
		cmap[i + 1] = (i + 1) / 2;
	}
	else
	{
		match[i] = i;
		cmap[i] = i / 2;
	}
  }
  cnvtxs = (nvtxs - 1) / 2 + 1;*/

  cudaDeviceSynchronize();
  gettimeofday(&end_part_cmatch,NULL);
  part_cmatch += (end_part_cmatch.tv_sec - begin_part_cmatch.tv_sec) * 1000 + (end_part_cmatch.tv_usec - begin_part_cmatch.tv_usec) / 1000.0;

  cudaDeviceSynchronize();
  gettimeofday(&begin_part_ccontract,NULL);
  hunyuangraph_cpu_create_cgraph(hunyuangraph_admin, graph, cnvtxs, match);
  cudaDeviceSynchronize();
  gettimeofday(&end_part_ccontract,NULL);
  part_ccontract += (end_part_ccontract.tv_sec - begin_part_ccontract.tv_sec) * 1000 + (end_part_ccontract.tv_usec - begin_part_ccontract.tv_usec) / 1000.0;

  return cnvtxs;
}

__global__ void init_vwgt(int *vwgt, int nvtxs)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;  

    if(ii < nvtxs)
        vwgt[ii] = 1;
}

__global__ void init_adjwgt(int *adjwgt, int nedges)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;  

    if(ii < nedges)
        adjwgt[ii] = 1;
}

/*Malloc and memcpy original graph from cpu to gpu*/
void hunyuangraph_malloc_original_coarseninfo(hunyuangraph_admin_t *hunyuangraph_admin,hunyuangraph_graph_t *graph)
{
    int nvtxs  = graph->nvtxs;
    int nedges = graph->nedges;
	graph->cuda_vwgt   = (int *)lmalloc_with_check(sizeof(int) * nvtxs,"vwgt");
	graph->cuda_xadj   = (int *)lmalloc_with_check(sizeof(int) * (nvtxs + 1),"xadj");
	graph->cuda_adjncy = (int *)lmalloc_with_check(sizeof(int) * nedges,"adjncy");
	graph->cuda_adjwgt = (int *)lmalloc_with_check(sizeof(int) * nedges,"adjwgt");

    // cudaMalloc((void**)&graph->cuda_xadj,(nvtxs+1)*sizeof(int));
    // cudaMalloc((void**)&graph->cuda_vwgt,nvtxs*sizeof(int));
    // cudaMalloc((void**)&graph->cuda_adjncy,nedges*sizeof(int));
    // cudaMalloc((void**)&graph->cuda_adjwgt,nedges*sizeof(int));

    cudaMemcpy(graph->cuda_xadj,graph->xadj,(nvtxs + 1) * sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(graph->cuda_adjncy,graph->adjncy,nedges*sizeof(int),cudaMemcpyHostToDevice);
    
    // ����CUDA��
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // ���õ�һ���˺���
    init_vwgt<<<(nvtxs + 127) / 128, 128, 0, stream>>>(graph->cuda_vwgt, nvtxs);
    // ���õڶ����˺���
    init_adjwgt<<<(nedges + 127) / 128, 128, 0, stream>>>(graph->cuda_adjwgt, nedges);

    // �ȴ����������??
    cudaStreamSynchronize(stream);

    // ����CUDA��
    cudaStreamDestroy(stream);

    // cudaMemcpy(graph->cuda_vwgt,graph->vwgt,nvtxs*sizeof(int),cudaMemcpyHostToDevice);
    // cudaMemcpy(graph->cuda_adjwgt,graph->adjwgt,nedges*sizeof(int),cudaMemcpyHostToDevice);
}

/*Malloc gpu coarsen graph params*/
void hunyuangraph_malloc_coarseninfo(hunyuangraph_admin_t *hunyuangraph_admin,hunyuangraph_graph_t *graph)
{
    int nvtxs  = graph->nvtxs;
    int nedges = graph->nedges;

    // cudaMalloc((void**)&graph->cuda_match,nvtxs * sizeof(int));
	graph->cuda_match = (int *)rmalloc_with_check(sizeof(int) * nvtxs,"match");
    // cudaMalloc((void**)&graph->cuda_cmap,nvtxs * sizeof(int));
	graph->cuda_cmap  = (int *)lmalloc_with_check(sizeof(int) * nvtxs,"cmap");
	graph->cuda_where = (int *)lmalloc_with_check(sizeof(int) * nvtxs,"where");		// 不可在k-way refinement再申请空间，会破坏栈的原则
}

/*CUDA-set each vertex pair adjacency list and weight params*/
__global__ void find_cnvtxsedge_original(int *txadj, int *xadj, int *match, int *adjncy,\
    int *cmap, int *tadjncy, int *tadjwgt, int *adjwgt, int nvtxs)
{
    // long long int i   = blockIdx.x * blockDim.x + threadIdx.x;
    int ii  = blockIdx.x * 4 + threadIdx.x / 32;
    int tid = threadIdx.x % 32;

    if(ii < nvtxs)
    {
        int u, t, pp, k, ptr, begin, end, iii;

        u = match[ii];
        t = cmap[ii];

        pp = txadj[t];
        if(ii > u)
        {
            begin = xadj[u];
            end   = xadj[u + 1];
            pp   += end - begin;
        }

        begin = xadj[ii];
        end   = xadj[ii + 1];

        for(iii = begin + tid, ptr = pp + tid;iii < end;iii += 32, ptr += 32)
        {
            k   = adjncy[iii];

            tadjncy[ptr] = cmap[k];
            tadjwgt[ptr] = adjwgt[iii];
        }
    }
}

__global__ void segment_sort(int *tadjncy, int *tadjwgt, int nedges, int *txadj, int cnvtxs)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < cnvtxs)
	{
		int begin, end, ptr, val;
		int i, j, k;
		begin = txadj[ii];
		end   = txadj[ii + 1];

		for(i = begin;i < end;i++)
		{
			ptr = i;
			val = tadjwgt[ptr];
			for(j = i + 1;j < end;j++)
				if(tadjwgt[j] < val) ptr = j, val = tadjwgt[ptr];
			val = tadjncy[ptr], tadjncy[ptr] = tadjncy[i], tadjncy[i] = val;
			val = tadjwgt[ptr], tadjwgt[ptr] = tadjwgt[i], tadjwgt[i] = val;
		}
	}
}

//Sort_cnedges2_part1<<<(cnvtxs + 3) / 4,128>>>(graph->tadjncy,graph->txadj,temp_scan,cnvtxs);
/*CUDA-Segmentation sorting part1-set scan array value 0 or 1*/
__global__ void Sort_cnedges2_part1(int *tadjncy, int *txadj, int *temp_scan, int cnvtxs)
{
    // long long int i   = blockIdx.x * blockDim.x + threadIdx.x;
	long long int ii  = blockIdx.x * 4 + threadIdx.x / 32;
	long long int tid = threadIdx.x % 32;

	if(ii < cnvtxs)
	{
		long long int j, begin, end, iii;

		begin  = txadj[ii];
		end    = txadj[ii + 1];

		for(iii = begin + tid;iii < end;iii += 32)
		{
			j   = tadjncy[iii];
				
			if(iii == begin)
			{
				if(j == ii) temp_scan[iii] = 0;
				else temp_scan[iii] = 1;
			}
			else 
			{
				if(j == ii) temp_scan[iii] = 0;
				else
				{
					if(j == tadjncy[iii - 1]) temp_scan[iii] = 0;
					else temp_scan[iii] = 1;
				}
			}
		}
	}
}

__global__ void Sort_cnedges2_part1_shared(int *tadjncy, int *txadj, int *temp_scan, int cnvtxs)
{
    int i   = blockIdx.x * blockDim.x + threadIdx.x;
    int ii  = i / 32;
    int tid = i - ii * 32;

    __shared__ int cache_txadj[8];
    // __shared__ int cache_tadjncy[4][32];
    int pid = ii % 4;

    if(ii < cnvtxs)
    {
        int j, begin, end;

        if(tid == 0) cache_txadj[pid] = txadj[ii];
        if(tid == 1) cache_txadj[pid + 1] = txadj[ii + 1];

        __syncthreads();

        // begin  = txadj[ii];
        // end    = txadj[ii + 1];

        //bank conflict
        begin = cache_txadj[pid];
        end   = cache_txadj[pid + 1];

        for(i = begin + tid;i < end;i += 32)
        {
            j   = tadjncy[i];
            // cache_tadjncy[pid][tid] = j;
            // __syncthreads();
            
            if(i == begin)
            {
                if(j == ii) temp_scan[i] = 0;
                else temp_scan[i] = 1;
            }
            else 
            {
                if(j == ii) temp_scan[i] = 0;
                else
                {
                    if(j == tadjncy[i - 1]) temp_scan[i] = 0;
                    // if(tid != 0 && j == cache_tadjncy[pid][tid - 1]) temp_scan[i] = 0;
                    // else if(tid == 0 && j == tadjncy[i - 1]) temp_scan[i] = 0;
                    else temp_scan[i] = 1;
                }
            }
        }
    }
}

/*CUDA-Segmentation sorting part2-set cxadj*/
__global__ void Sort_cnedges2_part2(int *txadj, int *temp_scan, int *cxadj, int cnvtxs)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

    if(ii < cnvtxs)
    { 
        int ppp = txadj[ii + 1];

        cxadj[ii + 1] = temp_scan[ppp - 1];
    }
    else if(ii == cnvtxs) cxadj[0] = 0;
} 

/*CUDA-Segmentation sorting part2.5-init cadjwgt and cadjncy*/
__global__ void Sort_cnedges2_part2_5(int *cadjwgt, int cnedges)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;  

    if(ii < cnedges)
        cadjwgt[ii] = 0;
}

/*CUDA-Segmentation sorting part3-deduplication and accumulation*/
__global__ void Sort_cnedges2_part3(int *tadjncy,int *txadj, int *tadjwgt,int *temp_scan, int *cxadj,int *cadjncy, int *cadjwgt, int cnvtxs)
{
    // long long int i   = blockIdx.x * blockDim.x + threadIdx.x;
    int ii  = blockIdx.x * 4 + threadIdx.x / 32;
    int tid = threadIdx.x % 32;

    // if(ii < cnvtxs)
    // {
    //     int ptr, j, begin, end;

    //     begin = txadj[ii];
    //     end   = txadj[ii + 1];

    //     for(i = begin + tid;i < end;i += 32)
    //     {
    //         j   = tadjncy[i];
    //         ptr = temp_scan[i] - 1;

    //         if(j != ii)
    //         {
    //             atomicAdd(&cadjwgt[ptr],tadjwgt[i]);

    //             if(i == begin) cadjncy[ptr] = j;
    //             else
    //             {
    //                 if(j != tadjncy[i - 1]) cadjncy[ptr] = j;
    //             }
    //         }
    //     }
    // }

	if(ii < cnvtxs)
	{
		int begin, end, j, k, iii;

		begin = txadj[ii];
        end   = txadj[ii + 1];

		for(iii = begin + tid;iii < end;iii += 32)
		{
			j = tadjncy[iii];
			k = temp_scan[iii] - 1;

			if(iii == begin)
			{
				if(j != ii)
				{
					cadjncy[k] = j;
					atomicAdd(&cadjwgt[k],tadjwgt[iii]);
				}
			}
			else 
			{
				if(j != ii)
				{
					if(j != tadjncy[iii - 1])
					{
						cadjncy[k] = j;
						atomicAdd(&cadjwgt[k],tadjwgt[iii]);
					}
					else atomicAdd(&cadjwgt[k],tadjwgt[iii]);
				}
			}
		}
	}
}

__global__ void exam_csr(int nvtxs, int *xadj, int *adjncy, int *adjwgt)
{
	for (int i = 0; i <= nvtxs; i++)
		printf("%d ", xadj[i]);

	printf("\nadjncy/adjwgt:\n");
	for (int i = 0; i < nvtxs; i++)
	{
		for (int j = xadj[i]; j < xadj[i + 1]; j++)
			printf("%d ", adjncy[j]);
		printf("\n");
		for (int j = xadj[i]; j < xadj[i + 1]; j++)
			printf("%d ", adjwgt[j]);
		printf("\n");
	}
	// printf("\n");
	// for (int i = 0; i < nvtxs; i++)
	// {
	// 	for (int j = xadj[i]; j < xadj[i + 1]; j++)
	// 		printf("%d ", adjwgt[j]);
	// 	printf("\n");
	// }
	// printf("\n");
}

/*Free cuda coarsen graph params*/
void hunyuangraph_free_coarsen(hunyuangraph_graph_t *graph)
{
  cudaFree(graph->cuda_match);
  cudaFree(graph->txadj);
  cudaFree(graph->tadjwgt);
  cudaFree(graph->tadjncy);
}

/*Create gpu coarsen graph by contract*/
void hunyuangraph_gpu_create_cgraph(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, hunyuangraph_graph_t *cgraph)
{
    int nvtxs  = graph->nvtxs;
    int nedges = graph->nedges;
    int cnvtxs = cgraph->nvtxs;

    cudaDeviceSynchronize();
    gettimeofday(&begin_exclusive_scan,NULL);
    thrust::exclusive_scan(thrust::device,graph->txadj,graph->txadj + cnvtxs + 1,graph->txadj);
    cudaDeviceSynchronize();
    gettimeofday(&end_exclusive_scan,NULL);
    sexclusive_scan += (end_exclusive_scan.tv_sec - begin_exclusive_scan.tv_sec) * 1000 + (end_exclusive_scan.tv_usec - begin_exclusive_scan.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_coarsen_malloc,NULL);
    // cudaMalloc((void**)&graph->tadjncy,nedges * sizeof(int));
    // cudaMalloc((void**)&graph->tadjwgt,nedges * sizeof(int));
	int *bb_keysB_d, *bb_valsB_d;
	bb_keysB_d = (int *)rmalloc_with_check(sizeof(int) * nedges,"bb_keysB_d");
	bb_valsB_d = (int *)rmalloc_with_check(sizeof(int) * nedges,"bb_valsB_d");
	graph->tadjncy = (int *)rmalloc_with_check(sizeof(int) * nedges,"tadjncy");
	graph->tadjwgt = (int *)rmalloc_with_check(sizeof(int) * nedges,"tadjwgt");
    cudaDeviceSynchronize();
    gettimeofday(&end_coarsen_malloc,NULL);
    coarsen_malloc += (end_coarsen_malloc.tv_sec - begin_coarsen_malloc.tv_sec) * 1000 + (end_coarsen_malloc.tv_usec - begin_coarsen_malloc.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_find_cnvtxsedge_original,NULL);
    find_cnvtxsedge_original<<<(nvtxs + 3) / 4,128>>>(graph->txadj,graph->cuda_xadj,graph->cuda_match,graph->cuda_adjncy,graph->cuda_cmap,\
        graph->tadjncy,graph->tadjwgt,graph->cuda_adjwgt,nvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&end_find_cnvtxsedge_original,NULL);
    sfind_cnvtxsedge_original += (end_find_cnvtxsedge_original.tv_sec - begin_find_cnvtxsedge_original.tv_sec) * 1000 + (end_find_cnvtxsedge_original.tv_usec - begin_find_cnvtxsedge_original.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_bb_segsort,NULL);
	int *bb_counter, *bb_id;
	// bb_keysB_d = (int *)rmalloc_with_check(sizeof(int) * nedges,"bb_keysB_d");
	// bb_valsB_d = (int *)rmalloc_with_check(sizeof(int) * nedges,"bb_valsB_d");
	bb_id      = (int *)rmalloc_with_check(sizeof(int) * cnvtxs,"bb_id");
	bb_counter = (int *)rmalloc_with_check(sizeof(int) * 13,"bb_counter");
    bb_segsort_me(graph->tadjncy, graph->tadjwgt, nedges, graph->txadj, cnvtxs, bb_counter, bb_id, bb_keysB_d, bb_valsB_d);
	// segment_sort<<<(cnvtxs + 127) / 128, 128>>>(graph->tadjncy, graph->tadjwgt, nedges, graph->txadj, cnvtxs);
	rfree_with_check(sizeof(int) * 13,"bb_counter");			//bb_counter
	rfree_with_check(sizeof(int) * cnvtxs,"bb_id");				//bb_id
	graph->tadjncy = bb_keysB_d;
	graph->tadjwgt = bb_valsB_d;
	rfree_with_check(sizeof(int) * nedges,"tadjwgt");			//tadjwgt
	rfree_with_check(sizeof(int) * nedges,"bb_keysB_d");		//tadjncy
    cudaDeviceSynchronize();
    gettimeofday(&end_bb_segsort,NULL);
    sbb_segsort += (end_bb_segsort.tv_sec - begin_bb_segsort.tv_sec) * 1000 + (end_bb_segsort.tv_usec - begin_bb_segsort.tv_usec) / 1000.0;

	// cudaDeviceSynchronize();
	// exam_csr<<<1, 1>>>(cnvtxs, graph->txadj, graph->tadjncy, graph->tadjwgt);
	// cudaDeviceSynchronize();

	// cudaDeviceSynchronize();
	// segment_sort<<<(cnvtxs + 127) / 128, 128>>>(graph->tadjncy, graph->tadjwgt, nedges, graph->txadj, cnvtxs);
	// cudaDeviceSynchronize();
	// cudaDeviceSynchronize();
	// exam_csr<<<1, 1>>>(cnvtxs, graph->txadj, graph->tadjncy, graph->tadjwgt);
	// cudaDeviceSynchronize();

    int *temp_scan;
    cudaDeviceSynchronize();
    gettimeofday(&begin_coarsen_malloc,NULL);
    // cudaMalloc((void**)&temp_scan, nedges * sizeof(int));
	temp_scan = (int *)rmalloc_with_check(sizeof(int) * nedges,"temp_scan");
    cudaDeviceSynchronize();
    gettimeofday(&end_coarsen_malloc,NULL);
    coarsen_malloc += (end_coarsen_malloc.tv_sec - begin_coarsen_malloc.tv_sec) * 1000 + (end_coarsen_malloc.tv_usec - begin_coarsen_malloc.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_Sort_cnedges_part1,NULL);
    Sort_cnedges2_part1<<<(cnvtxs + 3) / 4,128>>>(graph->tadjncy,graph->txadj,temp_scan,cnvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&end_Sort_cnedges_part1,NULL);
    sSort_cnedges_part1 += (end_Sort_cnedges_part1.tv_sec - begin_Sort_cnedges_part1.tv_sec) * 1000 + (end_Sort_cnedges_part1.tv_usec - begin_Sort_cnedges_part1.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_inclusive_scan,NULL);
	thrust::inclusive_scan(thrust::device,temp_scan,temp_scan + nedges,temp_scan);
    cudaDeviceSynchronize();
    gettimeofday(&end_inclusive_scan,NULL);
    sinclusive_scan2 += (end_inclusive_scan.tv_sec - begin_inclusive_scan.tv_sec) * 1000 + (end_inclusive_scan.tv_usec - begin_inclusive_scan.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_coarsen_malloc,NULL);
    // cudaMalloc((void**)&cgraph->cuda_xadj, (cnvtxs+1)*sizeof(int));
	cgraph->cuda_xadj = (int *)lmalloc_with_check(sizeof(int) * (cnvtxs + 1),"xadj");
    cudaDeviceSynchronize();
    gettimeofday(&end_coarsen_malloc,NULL);
    coarsen_malloc += (end_coarsen_malloc.tv_sec - begin_coarsen_malloc.tv_sec) * 1000 + (end_coarsen_malloc.tv_usec - begin_coarsen_malloc.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_Sort_cnedges_part2,NULL);
    Sort_cnedges2_part2<<<(cnvtxs + 128) / 128,128>>>(graph->txadj,temp_scan,cgraph->cuda_xadj,cnvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&end_Sort_cnedges_part2,NULL);
    sSort_cnedges_part2 += (end_Sort_cnedges_part2.tv_sec - begin_Sort_cnedges_part2.tv_sec) * 1000 + (end_Sort_cnedges_part2.tv_usec - begin_Sort_cnedges_part2.tv_usec) / 1000.0;

    cudaMemcpy(&cgraph->nedges, &cgraph->cuda_xadj[cnvtxs], sizeof(int), cudaMemcpyDeviceToHost); 

    cudaDeviceSynchronize();
    gettimeofday(&begin_coarsen_malloc,NULL);
    // cudaMalloc((void**)&cgraph->cuda_adjncy, cgraph->nedges * sizeof(int));
    // cudaMalloc((void**)&cgraph->cuda_adjwgt, cgraph->nedges * sizeof(int));
	cgraph->cuda_adjncy = (int *)lmalloc_with_check(sizeof(int) * cgraph->nedges,"adjncy");
	cgraph->cuda_adjwgt = (int *)lmalloc_with_check(sizeof(int) * cgraph->nedges,"adjwgt");
    cudaDeviceSynchronize();
    gettimeofday(&end_coarsen_malloc,NULL);
    coarsen_malloc += (end_coarsen_malloc.tv_sec - begin_coarsen_malloc.tv_sec) * 1000 + (end_coarsen_malloc.tv_usec - begin_coarsen_malloc.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_Sort_cnedges_part2_5,NULL);
    Sort_cnedges2_part2_5<<<(cgraph->nedges + 127) / 128,128>>>(cgraph->cuda_adjwgt,cgraph->nedges);
    cudaDeviceSynchronize();
    gettimeofday(&end_Sort_cnedges_part2_5,NULL);
    sSort_cnedges_part2_5 += (end_Sort_cnedges_part2_5.tv_sec - begin_Sort_cnedges_part2_5.tv_sec) * 1000 + (end_Sort_cnedges_part2_5.tv_usec - begin_Sort_cnedges_part2_5.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_Sort_cnedges_part3,NULL);
    Sort_cnedges2_part3<<<(cnvtxs + 3) / 4,128>>>(graph->tadjncy,graph->txadj,\
        graph->tadjwgt,temp_scan,cgraph->cuda_xadj,cgraph->cuda_adjncy,cgraph->cuda_adjwgt,cnvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&end_Sort_cnedges_part3,NULL);
    sSort_cnedges_part3 += (end_Sort_cnedges_part3.tv_sec - begin_Sort_cnedges_part3.tv_sec) * 1000 + (end_Sort_cnedges_part3.tv_usec - begin_Sort_cnedges_part3.tv_usec) / 1000.0;

	cgraph->tvwgt[0] = graph->tvwgt[0];  

	// printf("cnvtxs=%d cnedges=%d ",cnvtxs,cgraph->nedges);

	cudaDeviceSynchronize();
    gettimeofday(&begin_coarsen_free,NULL); 
	// rfree_with_check(sizeof(int) * nedges,"bb_valsB_d");		//bb_valsB_d
	// rfree_with_check(sizeof(int) * nedges,"bb_keysB_d");		//bb_keysB_d
	rfree_with_check(sizeof(int) * nedges,"temp_scan");			//temp_scan
	rfree_with_check(sizeof(int) * nedges,"tadjwgt");			//tadjwgt
	rfree_with_check(sizeof(int) * nedges,"tadjncy");			//tadjncy
	rfree_with_check(sizeof(int) * (cnvtxs + 1),"txadj");		//txadj
	rfree_with_check(sizeof(int) * nvtxs,"match");				//match
	// cudaFree(temp_scan);
    // cudaFree(graph->cuda_match);
    // cudaFree(graph->txadj);
    // cudaFree(graph->tadjwgt);
    // cudaFree(graph->tadjncy);
	cudaDeviceSynchronize();
    gettimeofday(&end_coarsen_free,NULL);
    coarsen_free += (end_coarsen_free.tv_sec - begin_coarsen_free.tv_sec) * 1000 + (end_coarsen_free.tv_usec - begin_coarsen_free.tv_usec) / 1000.0;
}

/*CUDA-init match array*/
__global__ void initcuda_match(int *match, int nvtxs)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

    if(ii < nvtxs)
        match[ii] = -1;
}

/*CUDA-hem matching*/
__global__ void cuda_hem(int nvtxs_hem, int *match, int *xadj, int *vwgt,int *adjwgt, int *adjncy, int maxvwgt_hem)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

    /*int maxvwgt, maxidx, maxwgt, i, j, ivwgt, k, jw;
    int ibegin, iend, begin, end;

    maxvwgt = maxvwgt_hem;

    if(ii < addition) ibegin = ii * (size + 1), iend = ibegin + size + 1;
    else ibegin = addition * (size + 1) + (ii - addition) * size, iend = ibegin + size;

    for(i = ibegin;i < iend;i++)
    {
        if(match[i] == -1)
        {
            begin = xadj[i];
            end   = xadj[i + 1];
            ivwgt = vwgt[i];

            maxidx = i;
            maxwgt = -1;

            for(j = begin;j < end;j++)
            {
                k  = adjncy[j];
                jw = adjwgt[j];
                if(match[k] == -1 && maxwgt < jw && ivwgt + vwgt[k] <= maxvwgt)
                {
                    maxidx = k;
                    maxwgt = jw;
                }
            }
                if(maxidx == i && 3 * ivwgt < maxvwgt)
                    maxidx = -1;

            if(maxidx != -1)
            {
                match[i] = maxidx;
                atomicExch(&match[maxidx],i);
            }
        }
    }*/

	int tt, nvtxs, maxvwgt, b, a, x, maxidx, maxwgt, i, j, ivwgt, k, jw;
  	int ibegin, iend, begin, end;

	tt      = 1024;
	nvtxs   = nvtxs_hem;
	maxvwgt = maxvwgt_hem;

	if(nvtxs % tt == 0)
	{
		b = nvtxs / tt;
		ibegin = ii * b;
		iend   = ibegin + b;
	}
	else 
	{
		b = nvtxs / tt;
		a = b + 1;
		x = nvtxs - b * tt;
		if(ii < x)
		{
			ibegin = ii * a;
			iend   = ibegin + a;
		}
		else
		{
			ibegin = ii * b + x;
			iend   = ibegin + b;
		}
	}
	for(i = ibegin;i < iend;i++)
	{
		if(match[i] == -1)
		{
			begin = xadj[i];
			end   = xadj[i + 1];
			ivwgt = vwgt[i];

			maxidx = i;
			maxwgt = -1;

			if (ivwgt < maxvwgt) 
			{
				for(j = begin;j < end;j++)
				{
					k  = adjncy[j];
					jw = adjwgt[j];
					if(match[k] == -1 && maxwgt < jw && ivwgt + vwgt[k] <= maxvwgt)
					{
						maxidx = k;
						maxwgt = jw;
					}
				}
				if(maxidx == i && 3 * ivwgt < maxvwgt)
					maxidx = -1;
			}

			if(maxidx != -1)
			{
				atomicCAS(&match[maxidx],-1,i);
				atomicExch(&match[i],maxidx);
			}
		}
	}
}

__global__ void cuda_hem1(int nvtxs, int *match, int *xadj, int *vwgt, int *adjwgt, int *adjncy, int maxvwgt)
{
  int pi;
  int ii;
  int i,j,k,maxidx,maxwgt;
  ii=blockIdx.x*blockDim.x+threadIdx.x;
  int b_start,b_end;
  int tt=1024;

  if(nvtxs%tt==0){
    b_start=ii*(nvtxs/tt);
    b_end=(ii+1)*(nvtxs/tt);
  }
  else{
    int b=nvtxs/tt;
    int a=b+1;
    int x=nvtxs-b*tt;

    if(ii<x){
      b_start=ii*a;
      b_end=(ii+1)*a;
    }
    else{
      b_start=x*a+(ii-x)*b;
      b_end=x*a+(ii+1-x)*b;
    }
  }

  for(pi=b_start;pi<b_end;pi++){
    i=pi;

    if(match[i]==-1){  
      maxidx=i;                                                                               
      maxwgt=-1;       

      for(j=xadj[i];j<xadj[i+1];j++){
        k=adjncy[j];

        if(match[k]==-1&&maxwgt<adjwgt[j]&&vwgt[i]+vwgt[k]<=maxvwgt){
          maxidx=k;
          maxwgt=adjwgt[j];
        }  
        if(maxidx==i&&3*vwgt[i]<maxvwgt){ 
          maxidx = -1;
        }
      }
      if(maxidx!=-1){    
        match[i] = maxidx;  
        atomicExch(&match[maxidx],i);                                 
      }
    }
  }
}

__global__ void cuda_hem_cpu(int nvtxs_hem, int *match, int *xadj, int *vwgt,int *adjwgt, int *adjncy, int maxvwgt_hem)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	int i,j,k,pi,last_unmatched,nvtxs,cnvtxs,maxidx,maxwgt,maxvwgt;
	int nunmatched = 0;

	nvtxs = nvtxs_hem;
	maxvwgt = maxvwgt_hem;

	for (cnvtxs=0, last_unmatched=0, pi=0; pi<nvtxs; pi++) 
	{
    	i = pi;

    if (match[i] == -1) {  /* Unmatched */
      maxidx = i;
      maxwgt = -1;

      if (vwgt[i] < maxvwgt) 
	  {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) { 
			if(pi >= last_unmatched) last_unmatched = pi + 1;
			else last_unmatched = last_unmatched + 1;
        //   last_unmatched = gk_max(pi, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            j = last_unmatched;
            if (match[j] == -1) 
			{
              maxidx = j;
              break;
            }
          }
        }
        else 
		{
          /* Find a heavy-edge matching, subject to maxvwgt constraints */
            /* single constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) 
			{
              k = adjncy[j];
              if (match[k] == -1 && maxwgt < adjwgt[j] && vwgt[i]+vwgt[k] <= maxvwgt) 
			  {
                maxidx = k;
                maxwgt = adjwgt[j];
              }
            }

            /* If it did not match, record for a 2-hop matching. */
            if (maxidx == i && 3*vwgt[i] < maxvwgt) {
              nunmatched++;
              maxidx = -1;
            }
        }
      }

      if (maxidx != -1) 
	  {
        // cmap[i]  = cmap[maxidx] = cnvtxs++;
        match[i] = maxidx;
        match[maxidx] = i;
      }
    }
  }
}

__global__ void cuda_hem_test(int nvtxs_hem, int *match, int *xadj, int *vwgt,int *adjwgt, int *adjncy, int maxvwgt_hem)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

    if(ii < nvtxs_hem)
	{
		if(ii % 2 == 0) match[ii] = ii + 1;
		else match[ii] = ii - 1;

		if(ii == nvtxs_hem - 1 && ii % 2 == 0) match[ii] = ii;
	}
}

__global__ void initcuda_mark(int *mark, int nvtxs)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

    if(ii < nvtxs)
        mark[ii] = 0;
}

__global__ void cuda_hem_229(int nvtxs, int *match, int *xadj, int *vwgt,int *adjwgt, int *adjncy, int maxvwgt, int *mark)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int maxidx, maxwgt, i, j, ivwgt, jw;
		int begin, end;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		ivwgt = vwgt[ii];

		int f = 0;
		if(match[ii] != -1 || ivwgt >= maxvwgt) ;
		else
		{
			maxidx = i;
			maxwgt = -1;

			for(i = begin;i < end;i++)
			{
				j = adjncy[i];
				if(j < i) f ++;
				else
				{
					jw = adjwgt[i];
					if(match[j] == -1 && maxwgt < jw && ivwgt + vwgt[j] <= maxvwgt)
					{
						maxidx = j;
						maxwgt = jw;
					}
				}
			}
			if(maxidx == i && 3 * ivwgt < maxvwgt)
				maxidx = -1;

			do
			{
				if(match[ii] != -1)
				{
					for(i = begin;i < end;i++)
					{
						j = adjncy[i];
						if(j > i) atomicAdd(&mark[j],1);
					}
				}
				else if(mark[ii] >= f) 
				{
					if(maxidx != -1 && atomicAdd(&match[maxidx],0) == -1)
					{
						atomicExch(&match[maxidx],ii);
						match[ii] = maxidx;
						// printf("ii=%d match=%d\n",ii,match[ii]);

						for(i = begin;i < end;i++)
						{
							j = adjncy[i];
							if(j > i) atomicAdd(&mark[j],1);
						}
					}
					else
					{
						if(maxidx != -1)
						{
							maxidx = i;
							maxwgt = -1;

							for(i = begin;i < end;i++)
							{
								j = adjncy[i];
								if(j < i) f ++;
								else
								{
									jw = adjwgt[i];
									if(match[j] == -1 && maxwgt < jw && ivwgt + vwgt[j] <= maxvwgt)
									{
										maxidx = j;
										maxwgt = jw;
									}
								}
							}
							if(maxidx == i && 3 * ivwgt < maxvwgt)
								maxidx = -1;
						}
					}
				}
			} while (mark[ii] < f && match[ii] == -1);
		}
	}
}

__global__ void cuda_hem_229_2(int nvtxs, int *match, int *xadj, int *vwgt,int *adjwgt, int *adjncy, int maxvwgt, int *mark, int *pp)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int maxidx, maxwgt, i, j, ivwgt, jw;
		int begin, end;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		ivwgt = vwgt[ii];

		int f = 0;
		// 是否已匹配
		if(match[ii] != -1) ;
		else
		{
			// 选择当前最佳匹配点
			maxidx = i;
			maxwgt = -1;

			for(i = begin;i < end;i++)
			{
				j = adjncy[i];
				if(j < i) f ++;
				else
				{
					jw = adjwgt[i];
					if(match[j] == -1 && maxwgt < jw && ivwgt + vwgt[j] <= maxvwgt)
					{
						maxidx = j;
						maxwgt = jw;
					}
				}
			}
			if(maxidx == i && 3 * ivwgt < maxvwgt) maxidx = -1;

			if(f == 0) pp[ii] = 1;

			// 当前点是否已解放
			if(atomicAdd(&mark[ii],0) >= f)
			{
				// 最佳匹配点是否存在
				if(maxidx != -1)
				{
					// 最佳匹配点是否已匹配
					atomicCAS(&match[maxidx],-1,ii);
					if(match[maxidx] == ii)
					{
						atomicExch(&match[ii],ii);
						// 解放邻点
						for(i = begin;i < end;i++)
						{
							j = adjncy[i];
							if(j > i) atomicAdd(&mark[j],1);
						}
						// 解放匹配点邻点
						begin = xadj[maxidx];
						end   = xadj[maxidx + 1];
						for(i = begin;i < end;i++)
						{
							j = adjncy[i];
							if(j > i) atomicAdd(&mark[j],1);
						}
					}
					// 最佳匹配点已匹配
					else
					{
						do
						{
							// 选择当前最佳匹配点
							maxidx = i;
							maxwgt = -1;

							for(i = begin;i < end;i++)
							{
								j = adjncy[i];
								if(j > i)
								{
									jw = adjwgt[i];
									if(match[j] == -1 && maxwgt < jw && ivwgt + vwgt[j] <= maxvwgt)
									{
										maxidx = j;
										maxwgt = jw;
									}
								}
							}
							if(maxidx == i && 3 * ivwgt < maxvwgt) maxidx = -1;

							// 最佳匹配点是否存在
							if(maxidx != -1)
							{
								// 最佳匹配点是否已匹配
								atomicCAS(&match[maxidx],-1,ii);
								if(match[maxidx] == ii) atomicExch(&match[ii],ii);
							}
							else
							{
								atomicExch(&match[ii],ii);
							}
						} while (match[maxidx] == -1);
						
						// 解放邻点
						for(i = begin;i < end;i++)
						{
							j = adjncy[i];
							if(j > i) atomicAdd(&mark[j],1);
						}
						// 若匹配成功则解放匹配点邻点
						if(maxidx != -1)
						{
							begin = xadj[maxidx];
							end   = xadj[maxidx + 1];
							for(i = begin;i < end;i++)
							{
								j = adjncy[i];
								if(j > i) atomicAdd(&mark[j],1);
							}
						}
					}
				}
				else //不存在，则匹配自身且解放邻点
				{
					atomicExch(&match[ii],ii);
					for(i = begin;i < end;i++)
					{
						j = adjncy[i];
						if(j > i) atomicAdd(&mark[j],1);
					}
				}

			}
		}
	}
}

__global__ void cuda_hem_229_3(int nvtxs, int *match, int *xadj, int *vwgt,int *adjwgt, int *adjncy, int maxvwgt)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int maxidx, maxwgt, j, k, ivwgt, jw;
		int begin, end;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		ivwgt = vwgt[ii];

		if(match[ii] == -1)
		{
			begin = xadj[ii];
			end   = xadj[ii + 1];
			ivwgt = vwgt[ii];

			maxidx = ii;
			maxwgt = -1;

			if (ivwgt < maxvwgt) 
			{
				for(j = begin;j < end;j++)
				{
					k  = adjncy[j];
					jw = adjwgt[j];
					if(match[k] == -1 && maxwgt < jw && ivwgt + vwgt[k] <= maxvwgt)
					{
						maxidx = k;
						maxwgt = jw;
					}
				}
				if(maxidx == ii && 3 * ivwgt < maxvwgt)
					maxidx = -1;
			}

			if(maxidx != -1)
			{
				atomicCAS(&match[maxidx],-1,ii);
				atomicExch(&match[ii],maxidx);
			}
		}
	}
}

__global__ void cuda_hem_301(int nvtxs, int *match, int *xadj, int *vwgt,int *adjwgt, int *adjncy, int maxvwgt, int *mark, int *pp)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int maxidx, maxwgt, j, k, ivwgt, jw;
		int begin, end;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		ivwgt = vwgt[ii];

		if(match[ii] == -1)
		{
			begin = xadj[ii];
			end   = xadj[ii + 1];
			// ivwgt = vwgt[ii];

			maxidx = ii;
			// maxwgt = -1;

			if (3 * ivwgt < maxvwgt) 
			{
				for(j = end - 1;j < begin;j--)
				{
					k  = adjncy[j];
					if(ivwgt + vwgt[k] <= maxvwgt)
					{
						if(atomicCAS(&match[k],-1,ii) && atomicCAS(&match[ii],-1,k))
							break;
					}
				}
			}
			else atomicExch(&match[ii],ii);
		}
	}
}

__global__ void reset_match(int nvtxs, int *match)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int t = atomicAdd(&match[ii],0);
		if(t != -1)
		{
			if(match[t] != ii) 
				atomicExch(&match[ii],-1);
		}
	}
}

/*CUDA-set conflict array*//*cuda_cleanv*/
/*CUDA-find cgraph vertex part1-remark the match array by s*//*findc1*/
/*CUDA-find cgraph vertex part2-make sure the pair small label vertex*//*findc2*/
__global__ void findc12_me(int *match, int *cmap, int nvtxs)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

    // if(ii < nvtxs)
    // {
    //     int t = match[ii];
    //     if((t != -1 && match[t] != ii) || t == -1)
    //     {
    //         match[ii] = ii;
    //         cmap[ii] = 1;
    //     }
    //     else 
    //     {
    //         if(ii <= t) cmap[ii] = 1;
    //         else cmap[ii] = 0;
    //     }

    //     if(ii == 0) cmap[0] = 0;
    // }
	if(ii < nvtxs)
	{
		int t = match[ii];
		if((t != -1 && match[t] != ii) || t == -1)
			match[ii] = ii;
		
		__syncthreads();

		if(ii == 0) cmap[ii] = 0;
		else 
		{
			t = match[ii];
			if(ii <= t) cmap[ii] = 1;
			else cmap[ii] = 0;
		}
	}
}

__global__ void findc1_me(int *match, int nvtxs)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if (ii < nvtxs)
	{
		int t = match[ii];
		if ((t != -1 && match[t] != ii) || t == -1)
			match[ii] = ii;
	}
}

__global__ void findc2_me(int *match, int *cmap, int nvtxs)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if (ii < nvtxs)
	{
		if (ii == 0)
			cmap[ii] = 0;
		else
		{
			int t = match[ii];
			if (ii <= t)
				cmap[ii] = 1;
			else
				cmap[ii] = 0;
		}
	}
}

/*CUDA-find cgraph vertex part4-make sure vertex pair real rdge*//*findc4*/
__global__ void findc4_cvwgt_me(int *match, int *cmap, int *txadj, int *xadj, int *cvwgt, int *vwgt, int nvtxs)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

    if(ii < nvtxs)
    {
        int u = match[ii];
        if(ii > u)
        {
            int t = cmap[u];
            cmap[ii] = t;
            cvwgt[t] = vwgt[ii] + vwgt[u];
        }
        else 
        {
            int t, begin, end, length; 
            t = cmap[ii];
            begin  = xadj[ii];
            end    = xadj[ii + 1];
            length = end - begin;
            if(u != ii)
            {
                begin = xadj[u];
                end   = xadj[u + 1];
                txadj[t] = length + end - begin;
            }
            else txadj[t] = length;

            if(ii == u) cvwgt[t] = vwgt[ii];
        }
    }
}

__global__ void exam_match(int nvtxs, int *match)
{
	for(int i = 0;i < nvtxs;i++)
		printf("i=%d match=%d\n",i,match[i]);
}

__global__ void exam_mark(int nvtxs, int *mark, int *pp)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

    if(ii < nvtxs)
	{
		pp[ii] = mark[ii];
	}
}

/*Get gpu graph matching params by hem*/
int hunyuangraph_gpu_match(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
    cudaDeviceSynchronize();
    gettimeofday(&begin_part_match,NULL);

    int nvtxs  = graph->nvtxs;
    int nedges = graph->nedges;
    int cnvtxs = 0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_initcuda_match,NULL);
    initcuda_match<<<(nvtxs + 127) / 128,128>>>(graph->cuda_match,nvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&end_initcuda_match,NULL);
    sinitcuda_match += (end_initcuda_match.tv_sec - begin_initcuda_match.tv_sec) * 1000 + (end_initcuda_match.tv_usec - begin_initcuda_match.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_cuda_match,NULL);
	// cuda_hem_229_3<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_match,graph->cuda_xadj,graph->cuda_vwgt,graph->cuda_adjwgt,graph->cuda_adjncy,\
	// 	hunyuangraph_admin->maxvwgt,mark,pp);
	// reset_match<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_match);
    // cuda_hem1<<<1024,1>>>(nvtxs,graph->cuda_match,graph->cuda_xadj,graph->cuda_vwgt,graph->cuda_adjwgt,graph->cuda_adjncy,\
        hunyuangraph_admin->maxvwgt);
	// cuda_hem_test<<<(nvtxs + 127) / 32,128>>>(nvtxs,graph->cuda_match,graph->cuda_xadj,graph->cuda_vwgt,graph->cuda_adjwgt,graph->cuda_adjncy,\
        hunyuangraph_admin->maxvwgt);
	// cuda_hem_cpu<<<1024,1>>>(nvtxs,graph->cuda_match,graph->cuda_xadj,graph->cuda_vwgt,graph->cuda_adjwgt,graph->cuda_adjncy,\
		hunyuangraph_admin->maxvwgt);
	for(int i = 0;i < 1;i++)
	{
		// cudaDeviceSynchronize();
		// cuda_hem_229_2<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_match,graph->cuda_xadj,graph->cuda_vwgt,graph->cuda_adjwgt,graph->cuda_adjncy,\
			hunyuangraph_admin->maxvwgt,mark,pp);
		cuda_hem_229_3<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_match,graph->cuda_xadj,graph->cuda_vwgt,graph->cuda_adjwgt,graph->cuda_adjncy,\
			hunyuangraph_admin->maxvwgt);
		// cudaDeviceSynchronize();
		reset_match<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_match);
		// cudaDeviceSynchronize();
		cuda_hem<<<1024,1>>>(nvtxs,graph->cuda_match,graph->cuda_xadj,graph->cuda_vwgt,graph->cuda_adjwgt,graph->cuda_adjncy,\
        	hunyuangraph_admin->maxvwgt);
	}
    cudaDeviceSynchronize();
    gettimeofday(&end_cuda_match,NULL);
    scuda_match += (end_cuda_match.tv_sec - begin_cuda_match.tv_sec) * 1000 + (end_cuda_match.tv_usec - begin_cuda_match.tv_usec) / 1000.0;

	// cudaDeviceSynchronize();
	// exam_match<<<1,1>>>(nvtxs,graph->cuda_match);
	// cudaDeviceSynchronize();

    cudaDeviceSynchronize();
    gettimeofday(&begin_findc2,NULL);
    // findc12_me<<<(nvtxs + 127) / 128,128>>>(graph->cuda_match,graph->cuda_cmap,nvtxs);
    findc1_me<<<(nvtxs + 127) / 128, 128>>>(graph->cuda_match, nvtxs);
    cudaDeviceSynchronize();
    findc2_me<<<(nvtxs + 127) / 128, 128>>>(graph->cuda_match, graph->cuda_cmap, nvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&end_findc2,NULL);
    sfindc2 += (end_findc2.tv_sec - begin_findc2.tv_sec) * 1000 + (end_findc2.tv_usec - begin_findc2.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_inclusive_scan2,NULL);
    thrust::inclusive_scan(thrust::device, graph->cuda_cmap, graph->cuda_cmap + nvtxs, graph->cuda_cmap);
    cudaDeviceSynchronize();
    gettimeofday(&end_inclusive_scan2,NULL);
    sinclusive_scan += (end_inclusive_scan2.tv_sec - begin_inclusive_scan2.tv_sec) * 1000 + (end_inclusive_scan2.tv_usec - begin_inclusive_scan2.tv_usec) / 1000.0;

    cudaMemcpy(&cnvtxs,&graph->cuda_cmap[nvtxs - 1], sizeof(int), cudaMemcpyDeviceToHost);
    cnvtxs++;

    hunyuangraph_graph_t *cgraph = hunyuangraph_set_gpu_cgraph(graph, cnvtxs); 
    cgraph->nvtxs = cnvtxs;

    cudaDeviceSynchronize();
    gettimeofday(&begin_coarsen_malloc,NULL);
    // cudaMalloc((void**)&graph->txadj, (cnvtxs + 1) * sizeof(int));
	graph->txadj = (int *)rmalloc_with_check(sizeof(int) * (cnvtxs + 1),"txadj");
    // cudaMalloc((void**)&cgraph->cuda_vwgt, cnvtxs * sizeof(int));
	cgraph->cuda_vwgt = (int *)lmalloc_with_check(sizeof(int) * cnvtxs,"vwgt");
    cudaDeviceSynchronize();
    gettimeofday(&end_coarsen_malloc,NULL);
    coarsen_malloc += (end_coarsen_malloc.tv_sec - begin_coarsen_malloc.tv_sec) * 1000 + (end_coarsen_malloc.tv_usec - begin_coarsen_malloc.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_findc4,NULL);
    findc4_cvwgt_me<<<(nvtxs + 127) / 128,128>>>(graph->cuda_match, graph->cuda_cmap, graph->txadj, graph->cuda_xadj, \
        cgraph->cuda_vwgt, graph->cuda_vwgt, nvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&end_findc4,NULL);
    sfindc4 += (end_findc4.tv_sec - begin_findc4.tv_sec) * 1000 + (end_findc4.tv_usec - begin_findc4.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&end_part_match,NULL);
    part_match += (end_part_match.tv_sec - begin_part_match.tv_sec) * 1000 + (end_part_match.tv_usec - begin_part_match.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_part_contract,NULL);
    hunyuangraph_gpu_create_cgraph(hunyuangraph_admin, graph, cgraph);
    cudaDeviceSynchronize();
    gettimeofday(&end_part_contract,NULL);
    part_contract += (end_part_contract.tv_sec - begin_part_contract.tv_sec) * 1000 + (end_part_contract.tv_usec - begin_part_contract.tv_usec) / 1000.0;
    
    return cnvtxs;
}

void hunyuangraph_memcpy_coarsentoinit(hunyuangraph_graph_t *graph)
{
    int nvtxs  = graph->nvtxs;
    int nedges = graph->nedges;

    graph->xadj   = (int *)malloc(sizeof(int) * (nvtxs + 1)); 
    graph->vwgt   = (int *)malloc(sizeof(int) * nvtxs); 
    graph->adjncy = (int *)malloc(sizeof(int) * nedges);
    graph->adjwgt = (int *)malloc(sizeof(int) * nedges);

	CPU_malloc(sizeof(int) * (nvtxs + 1));
	CPU_malloc(sizeof(int) * nvtxs);
	CPU_malloc(sizeof(int) * nedges);
	CPU_malloc(sizeof(int) * nedges);

    cudaMemcpy(graph->xadj, graph->cuda_xadj, (nvtxs + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(graph->vwgt, graph->cuda_vwgt, nvtxs * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(graph->adjncy, graph->cuda_adjncy, nedges * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(graph->adjwgt, graph->cuda_adjwgt, nedges * sizeof(int), cudaMemcpyDeviceToHost);
}

/*Gpu multilevel coarsen*/
hunyuangraph_graph_t *hunyuangarph_coarsen(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
    int level = 0;

    hunyuangraph_admin->maxvwgt = 1.5 * graph->tvwgt[0] / hunyuangraph_admin->Coarsen_threshold; 
    
    do
    {
		cudaDeviceSynchronize();
    	gettimeofday(&begin_coarsen_malloc,NULL);
        hunyuangraph_malloc_coarseninfo(hunyuangraph_admin,graph);
		cudaDeviceSynchronize();
		gettimeofday(&end_coarsen_malloc,NULL);
		coarsen_malloc += (end_coarsen_malloc.tv_sec - begin_coarsen_malloc.tv_sec) * 1000 + (end_coarsen_malloc.tv_usec - begin_coarsen_malloc.tv_usec) / 1000.0;
	
        hunyuangraph_gpu_match(hunyuangraph_admin,graph);

        graph = graph->coarser;

        level++;
        // printf("level=%d\n",level);
		// if(level == 1) break;
		// printf("state1=%d state2=%d state3=%d \n",graph->nvtxs > hunyuangraph_admin->Coarsen_threshold,\
			graph->nvtxs < 0.85 * graph->finer->nvtxs,\
			graph->nedges > graph->nvtxs / 2);
    }while(
		graph->nvtxs > hunyuangraph_admin->Coarsen_threshold && \
        graph->nvtxs < 0.85 * graph->finer->nvtxs && \
        graph->nedges > graph->nvtxs / 2); 

	cudaDeviceSynchronize();
    gettimeofday(&begin_coarsen_memcpy,NULL);
    hunyuangraph_memcpy_coarsentoinit(graph);  
	cudaDeviceSynchronize();
    gettimeofday(&end_coarsen_memcpy,NULL);
    coarsen_memcpy += (end_coarsen_memcpy.tv_sec - begin_coarsen_memcpy.tv_sec) * 1000 + (end_coarsen_memcpy.tv_usec - begin_coarsen_memcpy.tv_usec) / 1000.0;

    return graph;
}

/*Cpu multilevel coarsen*/
hunyuangraph_graph_t *hunyuangraph_cpu_coarsen(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
	int i, eqewgts, level=1;

	/* determine if the weights on the edges are all the same */
	for (eqewgts=1, i=1; i<graph->nedges; i++) {
		if (graph->adjwgt[0] != graph->adjwgt[i]) {
			eqewgts = 0;
			break;
		}
	}

	hunyuangraph_admin->maxvwgt = 1.5*graph->tvwgt[0]/hunyuangraph_admin->Coarsen_threshold;
 
	do{
		if(graph->cmap==NULL){
			graph->cmap=(int*)malloc(sizeof(int)*(graph->nvtxs));
			CPU_malloc(sizeof(int)*(graph->nvtxs));
		}

		// hunyuangraph_cpu_match_HEM(hunyuangraph_admin,graph,level);
		if (eqewgts || graph->nedges == 0)
          	hunyuangraph_cpu_match_RM(hunyuangraph_admin, graph);
        else
          	hunyuangraph_cpu_match_HEM(hunyuangraph_admin, graph);

		graph = graph->coarser;
		// eqewgts = 0;
		/*if(level == 1)
		{
			printf("xadj:\n");
			for(int i = 0;i <= graph->nvtxs;i++)
			{
				printf("%d ",graph->xadj[i]);
			}
			printf("\nvwgt:\n");
			for(int i = 0;i < graph->nvtxs;i++)
			{
				printf("%d ",graph->vwgt[i]);
			}
			printf("\nadjncy:\n");
			for(int i = 0;i < graph->nvtxs;i++)
			{
				for(int j = graph->xadj[i];j < graph->xadj[i + 1];j++)
					printf("%d ",graph->adjncy[j]);
				printf("\n");
			}
			printf("adjwgt:\n");
			for(int i = 0;i < graph->nvtxs;i++)
			{
				for(int j = graph->xadj[i];j < graph->xadj[i + 1];j++)
					printf("%d ",graph->adjwgt[j]);
				printf("\n");
			}
		}*/

		level++;

	}while(graph->nvtxs > hunyuangraph_admin->Coarsen_threshold && 
		graph->nvtxs < 0.85*graph->finer->nvtxs && 
		graph->nedges > graph->nvtxs/2);

  	return graph;
}

/*Malloc cpu 2way-refine params*/
void hunyuangraph_allocate_cpu_2waymem(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
  int nvtxs;
  nvtxs = graph->nvtxs;

  graph->pwgts=(int*)malloc(2*sizeof(int));
  graph->where=(int*)malloc(nvtxs*sizeof(int));
  graph->bndptr=(int*)malloc(nvtxs*sizeof(int));
  graph->bndlist=(int*)malloc(nvtxs*sizeof(int));
  graph->id=(int*)malloc(nvtxs*sizeof(int));
  graph->ed=(int*)malloc(nvtxs*sizeof(int));
}

/*Compute cpu 2way-refine params*/
void hunyuangraph_compute_cpu_2wayparam(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
  int i,j,nvtxs,nbnd,mincut,istart,iend,tid,ted,me;
  int *xadj,*vwgt,*adjncy,*adjwgt,*pwgts;
  int *where,*bndptr,*bndlist,*id,*ed;

  nvtxs= graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
  where=graph->where;
  id=graph->id;
  ed=graph->ed;

  pwgts=hunyuangraph_int_set_value(2,0,graph->pwgts);
  bndptr=hunyuangraph_int_set_value(nvtxs,-1,graph->bndptr);

  bndlist=graph->bndlist;

  for(i=0;i<nvtxs;i++){
    pwgts[where[i]] += vwgt[i];
  }

  for(nbnd=0,mincut=0,i=0;i<nvtxs;i++){
    istart=xadj[i];
    iend=xadj[i+1];
    me=where[i];
    tid=ted=0;

    for(j=istart;j<iend;j++){
      if(me==where[adjncy[j]]){
        tid+=adjwgt[j];
      }
      else{
        ted+=adjwgt[j];
      }
    }

    id[i]=tid;
    ed[i]=ted;

    if(ted>0||istart==iend){
      hunyuangraph_listinsert(nbnd,bndlist,bndptr,i);
      mincut+=ted;
    }
  }

  graph->mincut=mincut/2;
  graph->nbnd=nbnd; 

}

/*Compute cpu imbalance params*/
float hunyuangraph_compute_cpu_imbal(hunyuangraph_graph_t *graph, int nparts, float *part_balance, float *ubvec)
{
  int j,*pwgts;
  float max,cur;
  pwgts=graph->pwgts;
  max=-1.0;

  for(j=0;j<nparts;j++){
    cur=pwgts[j]*part_balance[j]-ubvec[0];

    if(cur>max){
      max=cur;
    }
  }

  return max;
}

/*Init queue */
 void hunyuangraph_queue_init(hunyuangraph_queue_t *queue, size_t maxnodes)
{
  int i;
  queue->nnodes=0;
  queue->maxnodes=maxnodes;
  queue->heap=(hunyuangraph_rkv_t*)malloc(sizeof(hunyuangraph_rkv_t)*maxnodes);
  queue->locator=(ssize_t*)malloc(sizeof(ssize_t)*maxnodes);
  CPU_malloc(sizeof(hunyuangraph_rkv_t) * maxnodes + sizeof(ssize_t) * maxnodes);

  for(i=0;i<maxnodes;i++){
    queue->locator[i]=-1;
  }

}

/*Create queue*/
hunyuangraph_queue_t *hunyuangraph_queue_create(size_t maxnodes)
{
  hunyuangraph_queue_t *queue; 
  queue = (hunyuangraph_queue_t *)malloc(sizeof(hunyuangraph_queue_t));
  CPU_malloc(sizeof(hunyuangraph_queue_t));

  hunyuangraph_queue_init(queue, maxnodes);

  return queue;
}

/*Insert node to queue*/
int hunyuangraph_queue_insert(hunyuangraph_queue_t *queue, int node, int key)
{
  ssize_t i,j;
  ssize_t *locator=queue->locator;
  hunyuangraph_rkv_t *heap=queue->heap;
  i = queue->nnodes++;

  while(i>0){
    j=(i-1)>>1;

    if(M_GT_N(key,heap[j].key)){
      heap[i]=heap[j];
      locator[heap[i].val]=i;
      i=j;
    }
    else
      break;
  }

  heap[i].key=key;
  heap[i].val=node;
  locator[node]=i;

  return 0;

}

/*Get top of queue*/
int hunyuangraph_queue_top(hunyuangraph_queue_t *queue)
{
  ssize_t i, j;
  ssize_t *locator;
  hunyuangraph_rkv_t *heap;

  int vtx, node;
  float key;

  if (queue->nnodes==0){
    return -1;
  }

  queue->nnodes--;
  heap=queue->heap;
  locator=queue->locator;
  vtx=heap[0].val;
  locator[vtx]=-1;

  if ((i=queue->nnodes)>0){
    key=heap[i].key;
    node=heap[i].val;
    i=0;

    while((j=2*i+1)<queue->nnodes){
      if(M_GT_N(heap[j].key,key)){
        if(j+1 < queue->nnodes&&M_GT_N(heap[j+1].key,heap[j].key)){
          j=j+1;
        }

        heap[i]=heap[j];
        locator[heap[i].val]=i;
        i=j;
      }
      else if(j+1<queue->nnodes&&M_GT_N(heap[j+1].key,key)){
        j=j+1;
        heap[i]=heap[j];
        locator[heap[i].val]=i;
        i=j;
      }
      else
        break;
    }

    heap[i].key=key;
    heap[i].val=node;
    locator[node]=i;

  }

  return vtx;

}

/*Delete node of queue*/
int hunyuangraph_queue_delete(hunyuangraph_queue_t *queue, int node)
{
  ssize_t i, j, nnodes;
  float newkey, oldkey;
  ssize_t *locator=queue->locator;

  hunyuangraph_rkv_t *heap=queue->heap;

  i=locator[node];
  locator[node]=-1;

  if(--queue->nnodes>0&&heap[queue->nnodes].val!=node) {
    node=heap[queue->nnodes].val;
    newkey=heap[queue->nnodes].key;
    oldkey=heap[i].key;

    if(M_GT_N(newkey,oldkey)){ 
      while(i>0){
        j=(i-1)>>1;

        if(M_GT_N(newkey,heap[j].key)){
          heap[i]=heap[j];
          locator[heap[i].val]=i;
          i=j;
        }
        else
          break;
      }
    }
    else{ 
      nnodes=queue->nnodes;

      while((j=(i<<1)+1)<nnodes){
        if(M_GT_N(heap[j].key,newkey)){
          if(j+1<nnodes&&M_GT_N(heap[j+1].key,heap[j].key)){
            j++;
          }

          heap[i]=heap[j];
          locator[heap[i].val]=i;
          i=j;
        }
        else if(j+1<nnodes&&M_GT_N(heap[j+1].key,newkey)){
          j++;
          heap[i]=heap[j];
          locator[heap[i].val]=i;
          i=j;
        }
        else
          break;
      }
    }

    heap[i].key=newkey;
    heap[i].val=node;
    locator[node]=i;

  }

  return 0;
}

/*Update queue node key*/
void hunyuangraph_queue_update(hunyuangraph_queue_t *queue, int node, int newkey)
{
  ssize_t i, j, nnodes;
  float oldkey;
  ssize_t *locator=queue->locator;

  hunyuangraph_rkv_t *heap=queue->heap;
  oldkey=heap[locator[node]].key;
  i=locator[node];

  if(M_GT_N(newkey,oldkey)){ 
    while(i>0){
      j=(i-1)>>1;

      if(M_GT_N(newkey,heap[j].key)){
        heap[i]=heap[j];
        locator[heap[i].val]=i;
        i=j;
      }
      else
        break;
    }
  }
  else{ 
    nnodes = queue->nnodes;

    while((j=(i<<1)+1)<nnodes){
      if(M_GT_N(heap[j].key,newkey)){
        if(j+1<nnodes&&M_GT_N(heap[j+1].key,heap[j].key)){
          j++;
        }

        heap[i]=heap[j];
        locator[heap[i].val]=i;
        i=j;
      }
      else if(j+1<nnodes&&M_GT_N(heap[j+1].key,newkey)){
        j++;
        heap[i]=heap[j];
        locator[heap[i].val]=i;
        i=j;
      }
      else
        break;
    }
  }

  heap[i].key=newkey;
  heap[i].val=node;
  locator[node]=i;
  return;

}

/*Free queue*/
void hunyuangraph_queue_free(hunyuangraph_queue_t *queue)
{
  if(queue == NULL) return;

  free(queue->heap);
  free(queue->locator);
  CPU_malloc(sizeof(queue->heap) + sizeof(queue->locator));

  queue->maxnodes = 0;

  free(queue);
  CPU_malloc(sizeof(hunyuangraph_queue_t));
}

/*Reset queue*/
void hunyuangraph_queue_reset(hunyuangraph_queue_t *queue)
{
  ssize_t i;
  ssize_t *locator=queue->locator;

  hunyuangraph_rkv_t *heap=queue->heap;

  for(i=queue->nnodes-1;i>=0;i--){
    locator[heap[i].val]=-1;
  }

  queue->nnodes=0;

}

/*Balance two partition by moving boundary vertex*/
void hunyuangraph_bndvertex_2way_bal(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, float *ntpwgts)
{
  int i,ii,j,k,kwgt,nvtxs,nbnd,nswaps,from,to,temp;
  int *xadj,*vwgt,*adjncy,*adjwgt,*where,*id,*ed,*bndptr,*bndlist,*pwgts;
  int *moved,*perm;

  hunyuangraph_queue_t *queue;
  int higain,mincut,mindiff;
  int tpwgts[2];

  nvtxs=graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
  where=graph->where;
  id=graph->id;
  ed=graph->ed;
  pwgts=graph->pwgts;
  bndptr=graph->bndptr;
  bndlist=graph->bndlist;

  moved=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  perm=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  CPU_malloc(sizeof(int) * nvtxs);

  tpwgts[0]=graph->tvwgt[0]*ntpwgts[0];
  tpwgts[1]=graph->tvwgt[0]-tpwgts[0];
  mindiff=abs(tpwgts[0]-pwgts[0]);
  from=(pwgts[0]<tpwgts[0]?1:0);
  to=(from+1)%2;

  queue=hunyuangraph_queue_create(nvtxs);
  hunyuangraph_int_set_value(nvtxs,-1,moved);
  nbnd=graph->nbnd;
  hunyuangraph_int_randarrayofp(nbnd,perm,nbnd/5,1);

  for(ii=0;ii<nbnd;ii++){
    i=perm[ii];

    if(where[bndlist[i]]==from&&vwgt[bndlist[i]]<=mindiff){
      hunyuangraph_queue_insert(queue,bndlist[i],ed[bndlist[i]]-id[bndlist[i]]);
    }
  }

  mincut=graph->mincut;

  for(nswaps=0;nswaps<nvtxs;nswaps++) 
  {
    if((higain=hunyuangraph_queue_top(queue))==-1)
      break;
    if(pwgts[to]+vwgt[higain]>tpwgts[to])
      break;

    mincut-=(ed[higain]-id[higain]);
    hunyuangraph_add_sub(pwgts[to],pwgts[from],vwgt[higain]);

    where[higain]=to;
    moved[higain]=nswaps;
    hunyuangraph_swap(id[higain],ed[higain],temp);

    if(ed[higain]==0&&xadj[higain]<xadj[higain+1]){ 
      hunyuangraph_listdelete(nbnd,bndlist,bndptr,higain);
    }

    for(j=xadj[higain];j<xadj[higain+1];j++){
      k=adjncy[j];
      kwgt=(to==where[k]?adjwgt[j]:-adjwgt[j]);
      hunyuangraph_add_sub(id[k],ed[k],kwgt);

      if(bndptr[k]!=-1){ 
        if(ed[k]==0){ 
          hunyuangraph_listdelete(nbnd,bndlist,bndptr,k);

          if(moved[k]==-1&&where[k]==from&&vwgt[k]<=mindiff){ 
            hunyuangraph_queue_delete(queue,k);
          }
        }
        else{ 
          if(moved[k]==-1&&where[k]==from&&vwgt[k]<=mindiff){
            hunyuangraph_queue_update(queue,k,ed[k]-id[k]);
          }
        }
      }
      else{
        if(ed[k]>0){  
          hunyuangraph_listinsert(nbnd,bndlist,bndptr,k);

          if(moved[k]==-1&&where[k]==from&&vwgt[k]<=mindiff){ 
            hunyuangraph_queue_insert(queue,k,ed[k]-id[k]);
          }
        }
      }
    }
  }

  graph->mincut=mincut;
  graph->nbnd=nbnd;
  hunyuangraph_queue_free(queue);

}

/*Balance 2-way partition*/
void hunyuangraph_2way_bal(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, float *ntpwgts)
{
  if(hunyuangraph_compute_cpu_imbal(graph,2,hunyuangraph_admin->part_balance,hunyuangraph_admin->ubfactors)<=0){ 
    return;
  }

  if(abs(ntpwgts[0]*graph->tvwgt[0]-graph->pwgts[0])<3*graph->tvwgt[0]/graph->nvtxs){
    return;
  }

  hunyuangraph_bndvertex_2way_bal(hunyuangraph_admin,graph,ntpwgts);
}

__global__ void init_moveto(int *where, int *moveto, int nvtxs)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nvtxs)
  {
    int t = where[ii];
    moveto[ii] = t;
  }
}

__global__ void update_moveto(int nvtxs, int nedges, int *xadj, int *adjncy, int *adjwgt, int *where, int *moveto, int *list, int *gainv)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nvtxs)
  {
    int i, j, ewj, wj, id, ed, g, moveto_j;
    int me, other, begin, end;
    me    = where[ii];
    other = (me + 1) % 2;
    begin = xadj[ii];
    end   = xadj[ii + 1];
    id    = 0;
    ed    = 0; 
    
    for(i = begin;i < end;i++)
    {
      j   = adjncy[i];
      ewj = adjwgt[i];
      wj  = where[j];

      if(me == wj) id += ewj;
      else ed += ewj;
    }

    g = ed - id;
    gainv[ii] = g;

    // 怎么筛，筛完怎么记录，要不要记录gainv
      // 移动�?1，不移动�?0
    // if(g > -0.25 * id) printf("ii1=%d\n",ii), list[ii] = 1, moveto[ii] = other;
    if(g > -0.25 * id) list[ii] = 1, moveto[ii] = other;
    else list[ii] = 0;
    // printf("ii1=%d moveto=%d\n",ii,moveto[ii]);

    __syncthreads();  //什么同步？？？！！！！

    // 晒第二遍 怎么�?
    if(g > -0.25 * id)
    {
      g = 0;
      for(i = begin;i < end;i++)
      {
        j   = adjncy[i];
        ewj = adjwgt[i];
        wj  = where[j];
        moveto_j = moveto[j];

        if(list[j] == 1 && (gainv[j] > gainv[ii] || gainv[j] == gainv[ii] && j < ii))
          moveto_j = wj;
        if(moveto_j == other) g += ewj;
        else g -= ewj;
      }

      // 需要吗�?
      // gainv[ii] = g;

      // if(g > -0.25 * id) printf("ii2yes=%d\n",ii);
      // else printf("ii2=%d me=%d to=%d\n",ii,me,other), list[ii] = 0, moveto[ii] = me;
      if(g > -0.25 * id) ;
      else list[ii] = 0, moveto[ii] = me;
    }

    __syncthreads();

    // 移动顶点，上面的操作对每层图要进行一次还是多次，什么时候结束，怎么判断 
      //先只筛一�?
    if(list[ii] == 1) where[ii] = moveto[ii];

    // __syncthreads();

    // 如何使得分区平衡 

    /*for(i = begin;i < end;i++)
    {
      j   = adjncy[i];
      ewj = adjwgt[i];
      wj  = where[j];

      if(me == wj) id += ewj;
      else ed += ewj;
    }*/ // 已存放到gainv�?

  }
}

__global__ void filter_first(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *moveto, int *list, int *gainv)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, ewj, wj, id, ed, g;
		int me, other, begin, end;
		me    = where[ii];
		other = (me + 1) % 2;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		id    = 0;
		ed    = 0; 
		
		for(i = begin;i < end;i++)
		{
			j   = adjncy[i];
			ewj = adjwgt[i];
			wj  = where[j];

			if(me == wj) id += ewj;
			else ed += ewj;
		}

		g = ed - id;
		gainv[ii] = g;

		// 怎么筛，筛完怎么记录，要不要记录gainv
		// 移动�?1，不移动�?0
		if(g > -0.25 * id) list[ii] = 1, moveto[ii] = other;
		else list[ii] = 0;

	}
}

__global__ void filter_second(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *moveto, int *list, int *gainv)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, ewj, wj, l, g0, g, gj, moveto_j;
		int me, other, begin, end;
		me    = where[ii];
		other = (me + 1) % 2;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		l     = list[ii];
		g0    = gainv[ii];

		// 晒第二遍 怎么�?
		if(l)
		{
			g = 0;
			for(i = begin;i < end;i++)
			{
				j   = adjncy[i];
				ewj = adjwgt[i];
				wj  = where[j];
				gj  = gainv[j];
				moveto_j = moveto[j];

				if(list[j] == 1 && (gj > g0 || gj == g0 && j < ii))
					moveto_j = wj;
				if(moveto_j == other) g += ewj;
				else g -= ewj;
			}

			if(g > 0) where[ii] = other;
			else list[ii] = 0, moveto[ii] = me;
		}
	}
}

__global__ void filter_first_atomic(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *moveto, int *list, int *gainv)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, ewj, wj, id, ed, g;
		int me, other, begin, end;
		me    = where[ii];
		other = (me + 1) % 2;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		id    = 0;
		ed    = 0; 
		
		for(i = begin;i < end;i++)
		{
			j   = adjncy[i];
			ewj = adjwgt[i];
			wj  = where[j];

			if(me == wj) id += ewj;
			else ed += ewj;
		}

		g = ed - id;
		gainv[ii] = g;

		// 怎么筛，筛完怎么记录，要不要记录gainv
		// 移动�?1，不移动�?0
		if(g > -0.25 * id) list[ii] = 1, moveto[ii] = other;
		else list[ii] = 0;

	}
}

__global__ void filter_second_atomic(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *moveto, int *list, int *gainv)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, ewj, wj, l, g0, g, gj, moveto_j;
		int me, other, begin, end;
		me    = where[ii];
		other = (me + 1) % 2;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		l     = list[ii];
		g0    = gainv[ii];

		// 晒第二遍 怎么�?
		if(l)
		{
			g = 0;
			for(i = begin;i < end;i++)
			{
				j   = adjncy[i];
				ewj = adjwgt[i];
				wj  = where[j];
				gj  = gainv[j];
				moveto_j = moveto[j];

				if(list[j] == 1 && (gj > g0 || gj == g0 && j < ii))
					moveto_j = wj;
				if(moveto_j == other) g += ewj;
				else g -= ewj;
			}

			if(g > 0) where[ii] = other;
			else list[ii] = 0, moveto[ii] = me;
		}
	}
}

__global__ void filter_atomic(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, me, other, k, wk, ptr, begin, end, ed, id;
		// for(i = ii;i < nvtxs;i += 1024)
		// {
			ptr = ii;
			// if(ptr >= nvtxs) break;
				
			begin = xadj[ptr];
			end   = xadj[ptr + 1];
			me    = where[ptr];
			other = (me + 1) % 2;
			ed = 0;
			id = 0;

			for(j = begin;j < end;j++)
			{
				k  = adjncy[j];
				wk = where[k];

				if(wk == me) ed += adjwgt[j];
				else id += adjwgt[j];
			}

			// printf("me=%d other=%d\n",me,other);
			// if(me + other == 1 && me != other) ;
			// else printf("No\n");
			// if(ed > id) atomicExch(&where[ptr], other);
			if(ed <= id) 
			{
				int oldVal = atomicExch(&where[ptr], other);
			}
			// if(ed > id) 
			// {
			// 	if(me == 0) where[ptr] = 1;
			// 	else where[ptr] = 0;
			// }
		// }
	}
}

__global__ void filter_atomic_group(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where)
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;

	// if(ii < nvtxs)
	// {
		int i, j, me, other, k, wk, ptr, begin, end, ed, id;
		for(i = ii;i < nvtxs;i += 8)
		{
			ptr = i;
			if(ptr >= nvtxs) break;
				
			begin = xadj[ptr];
			end   = xadj[ptr + 1];
			me    = where[ptr];
			other = (me + 1) % 2;
			ed = 0;
			id = 0;

			for(j = begin;j < end;j++)
			{
				k  = adjncy[j];
				wk = where[k];

				if(wk == me) ed += adjwgt[j];
				else id += adjwgt[j];
			}

			// if(me + other == 1 && me != other) ;
			// else printf("No\n");
			// if(ed > id) atomicExch(&where[ptr], me);
			if(ed <= id) 
			{
				int oldVal = atomicExch(&where[ptr], other);
			}

		}
	// }
}

__global__ void update_moveto_atomic(int nvtxs, int nedges, int *xadj, int *adjncy, int *adjwgt, int *where, int *moveto, int *list, int *gainv)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nvtxs)
  {
    int i, j, ewj, wj, id, ed, g, moveto_j;
    int me, other, begin, end;
    me    = where[ii];
    other = (me + 1) % 2;
    begin = xadj[ii];
    end   = xadj[ii + 1];
    id    = 0;
    ed    = 0; 
    
    for(i = begin;i < end;i++)
    {
      j   = adjncy[i];
      ewj = adjwgt[i];
      wj  = where[j];

      if(me == wj) id += ewj;
      else ed += ewj;
    }

    g = ed - id;
    gainv[ii] = g;

    // 怎么筛，筛完怎么记录，要不要记录gainv
      // 移动�?1，不移动�?0
    // if(g > -0.25 * id) printf("ii1=%d\n",ii), list[ii] = 1, moveto[ii] = other;
    if(g > -0.25 * id) list[ii] = 1, moveto[ii] = other;
    else list[ii] = 0;
    // printf("ii1=%d moveto=%d\n",ii,moveto[ii]);

    __syncthreads();

    // 晒第二遍 怎么�?
    if(g > -0.25 * id)
    {
      g = 0;
      for(i = begin;i < end;i++)
      {
        j   = adjncy[i];
        ewj = adjwgt[i];
        wj  = where[j];
        moveto_j = moveto[j];

        if(list[j] > 1 && (gainv[j] > gainv[ii] || gainv[j] == gainv[ii] && j < ii))
          moveto_j = wj;
        if(moveto_j == other) g += ewj;
        else g -= ewj;
      }

      // 需要吗�?
      // gainv[ii] = g;

      // if(g > -0.25 * id) printf("ii2yes=%d\n",ii);
      // else printf("ii2=%d me=%d to=%d\n",ii,me,other), list[ii] = 0, moveto[ii] = me;
      if(g > -0.25 * id) ;
      else list[ii] = 0, moveto[ii] = me;
    }

    __syncthreads();

    // 移动顶点，上面的操作对每层图要进行一次还是多次，什么时候结束，怎么判断 
      //先只筛一�?
    if(list[ii] == 1) where[ii] = moveto[ii];

  }
}

__global__ void compute_pwgts(int nvtxs, int *vwgt, int *where, int *gainv)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 0) gainv[ii] = vwgt[ii];
		else gainv[ii] = 0;
	}
}

__device__ void warpReduction6(volatile int *share_num, int tid ,int blocksize)
{
    if(blocksize >= 64) share_num[tid] += share_num[tid + 32];
    if(blocksize >= 32) share_num[tid] += share_num[tid + 16];
    if(blocksize >= 16) share_num[tid] += share_num[tid + 8];
    if(blocksize >= 8) share_num[tid] += share_num[tid + 4];
    if(blocksize >= 4) share_num[tid] += share_num[tid + 2];
    if(blocksize >= 2) share_num[tid] += share_num[tid + 1];
}

__global__ void reduction6(int *num, int length)
{
    int ii  = blockIdx.x * (blockDim.x * 2) + threadIdx.x;
    int tid = threadIdx.x;
    
    extern __shared__ int share_num[];

    if(ii + blockDim.x < length) share_num[tid] = num[ii] + num[ii + blockDim.x];
    else if(ii < length) share_num[tid] = num[ii];
    else share_num[tid] = 0;

     __syncthreads();

    if(ii < length)
    {
        if(blockDim.x >= 512) 
        {
            if(tid < 256) share_num[tid] += share_num[tid + 256];
            __syncthreads();
        }
        if(blockDim.x >= 256) 
        {
            if(tid < 128) share_num[tid] += share_num[tid + 128];
            __syncthreads();
        }
        if(blockDim.x >= 128) 
        {
            if(tid < 64) share_num[tid] += share_num[tid + 64];
            __syncthreads();
        }

        if(tid < 32) warpReduction6(share_num, tid, blockDim.x);

        if(tid == 0) num[blockIdx.x] = share_num[0];
    }
}

__global__ void compute_gainkp(int nvtxs, int nedges, int from, int *xadj, int *adjncy, int *adjwgt, int *where, kp_t *kp)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, ewj, wj, id, ed, g;
		int me, begin, end;
		me    = where[ii];
		begin = xadj[ii];
		end   = xadj[ii + 1];
		id    = 0;
		ed    = 0; 
		
		for(i = begin;i < end;i++)
		{
			j   = adjncy[i];
			ewj = adjwgt[i];
			wj  = where[j];

			if(me == wj) id += ewj;
			else ed += ewj;
		}

		g = ed - id;
		if(from == me) kp[ii].key = g;
		else kp[ii].key = -nedges;
		kp[ii].ptr = ii;
	}
}

struct compRule
{
	__host__ __device__ bool operator()(const kp_t &p1,const kp_t &p2)
	{
		if (p1.key != p2.key)
			return p1.key > p2.key;
		else
      return p1.ptr < p2.ptr;
	}
};

__global__ void compute_gainvkp(int nvtxs, int *vwgt, kp_t *kp, int *gainv)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nvtxs)
  {
    int ptr = kp[ii].ptr;
    gainv[ii] = vwgt[ptr];
  }
}

__global__ void rebalancekp(int nvtxs, int num, int to, int *where, kp_t *kp, int *gainv)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nvtxs)
  {
    int ptr, prefixsum;
    ptr = kp[ii].ptr;
    prefixsum = gainv[ii];
    if(prefixsum <= num) where[ptr] = to;
  }
}

__global__ void exam_csr(int nvtxs, int *xadj, int *adjncy)
{
	for (int i = 0; i <= nvtxs; i++)
		printf("%d ", xadj[i]);

	printf("\nadjncy:\n");
	for (int i = 0; i < nvtxs; i++)
	{
		for (int j = xadj[i]; j < xadj[i + 1]; j++)
			printf("%d ", adjncy[j]);
		printf("\n");
	}
}

__global__ void exam_where(int nvtxs, int *where)
{
	for(int i = 0;i < nvtxs;i++)
		printf("%d ",where[i]);
	printf("\n");
}

__global__ void exam_gain_gainptr(int nvtxs, int *gain, int *gain_ptr)
{
  printf("gain:");
  for(int i = 0;i < nvtxs;i++)
    printf("%d ",gain[i]);
  printf("\n");
  printf("gain_ptr:");
  for(int i = 0;i < nvtxs;i++)
    printf("%d ",gain_ptr[i]);
  printf("\n");
}

__global__ void exam_kp(int nvtxs, kp_t *kp)
{
  printf("gain:");
  for(int i = 0;i < nvtxs;i++)
    printf("%d ",kp[i].key);
  printf("\n");
  printf("gain_ptr:");
  for(int i = 0;i < nvtxs;i++)
    printf("%d ",kp[i].ptr);
  printf("\n");
}

int compute_edgecut(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where)
{
  int edgecut, me;
  edgecut = 0;
  for(int i = 0;i < nvtxs;i++)
  {
    me =  where[i];
    for(int j = xadj[i];j < xadj[i + 1];j++)
      if(where[adjncy[j]] != me) edgecut += adjwgt[j];
  }
  return edgecut / 2;
}

__global__ void calculate_edgecut(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *temp)
{
	/*int iii = threadIdx.x;
	int ii  = blockIdx.x * 4 + iii / 32;
	int tid = iii % 32;

	extern __shared__ int cache_d[128];
	cache_d[iii] = 0;

	if (ii < nvtxs)
	{
		int begin, end, me, i, j;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		me    = where[ii];

		for(i = begin + tid;i < end;i += 32)
		{
			j = adjncy[i];
			if(where[j] != me) cache_d[iii] += adjwgt[i];
		}
	}
	__syncthreads();

	// Block内求和
	if(iii < 64) cache_d[iii] += cache_d[iii + 64];
	__syncthreads();
	if(iii < 32) cache_d[iii] += cache_d[iii + 32];
	__syncthreads();
	// Warp内求和
	if(iii < 16) cache_d[iii] += cache_d[iii + 16];
	if(iii < 8) cache_d[iii] += cache_d[iii + 8];
	if(iii < 4) cache_d[iii] += cache_d[iii + 4];
	if(iii < 2) cache_d[iii] += cache_d[iii + 2];
	if(iii < 1) cache_d[iii] += cache_d[iii + 1];
	
	if(iii == 0) temp[blockIdx.x] = cache_d[0];*/

	int ii  = blockIdx.x;
	int iii = threadIdx.x;

	__shared__ int cache_d[32];
	cache_d[iii] = 0;

	int begin, end, me, i, j;
	begin = xadj[ii];
	end   = xadj[ii + 1];
	me    = where[ii];

	for(i = begin + iii;i < end;i += 32)
	{
		j = adjncy[i];
		if(where[j] != me) cache_d[iii] += adjwgt[i];
	}
	__syncthreads();

	// Warp内求和
	if(iii < 16) cache_d[iii] += cache_d[iii + 16];
	__syncthreads();
	if(iii < 8) cache_d[iii] += cache_d[iii + 8];
	__syncthreads();
	if(iii < 4) cache_d[iii] += cache_d[iii + 4];
	__syncthreads();
	if(iii < 2) cache_d[iii] += cache_d[iii + 2];
	__syncthreads();
	if(iii < 1) cache_d[iii] += cache_d[iii + 1];
	__syncthreads();
	
	if(iii == 0) temp[ii] = cache_d[0];
}

__global__ void heat(int nvtxs, int *h)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if (ii < nvtxs)
	{
		int t = h[ii];
	}
}

void jetLP(hunyuangraph_graph_t *graph)
{
  int nvtxs, from, to;
  int *move_to, *list, *gainv;
  kp_t *kp;

  nvtxs = graph->nvtxs;

  cudaDeviceSynchronize();
  gettimeofday(&begin_malloc_2way, NULL);
  cudaMalloc((void**)&move_to,sizeof(int) * nvtxs);
  cudaMalloc((void**)&list,sizeof(int) * nvtxs);
  cudaMalloc((void**)&gainv,sizeof(int) * nvtxs);
  cudaMalloc((void**)&kp,sizeof(kp_t) * nvtxs);
  cudaDeviceSynchronize();
  gettimeofday(&end_malloc_2way, NULL);
//   malloc_2way += (end_malloc_2way.tv_sec - begin_malloc_2way.tv_sec) * 1000 + (end_malloc_2way.tv_usec - begin_malloc_2way.tv_usec) / 1000.0;

  cudaDeviceSynchronize();
  gettimeofday(&begin_initmoveto, NULL);
  init_moveto<<<(nvtxs + 127) / 128,128>>>(graph->cuda_where,move_to,nvtxs);
  cudaDeviceSynchronize();
  gettimeofday(&end_initmoveto, NULL);
  initmoveto += (end_initmoveto.tv_sec - begin_initmoveto.tv_sec) * 1000 + (end_initmoveto.tv_usec - begin_initmoveto.tv_usec) / 1000.0;

  cudaDeviceSynchronize();
  gettimeofday(&begin_updatemoveto, NULL);
  for(int i = 0;i < 1;i++)
  {
  	// filter_first<<<(nvtxs + 127) / 128,128>>>(nvtxs, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where, move_to, list, gainv);
  	// filter_second<<<(nvtxs + 127) / 128,128>>>(nvtxs, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where, move_to, list, gainv);
	
	// filter_first_atomic<<<(nvtxs + 127) / 128,128>>>(nvtxs, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where, move_to, list, gainv);
  	// filter_second_atomic<<<(nvtxs + 127) / 128,128>>>(nvtxs, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where, move_to, list, gainv);
	filter_atomic<<<(nvtxs + 127) / 128,128>>>(nvtxs, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where);
	// filter_atomic_group<<<(nvtxs / 8) / 32 + 1,32>>>(nvtxs, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where);

	// update_moveto_atomic<<<(nvtxs + 127) / 128,128>>>(nvtxs, graph->nedges, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where, move_to, list, gainv);
  }
  cudaDeviceSynchronize();
  gettimeofday(&end_updatemoveto, NULL);
  updatemoveto += (end_updatemoveto.tv_sec - begin_updatemoveto.tv_sec) * 1000 + (end_updatemoveto.tv_usec - begin_updatemoveto.tv_usec) / 1000.0;

  cudaDeviceSynchronize();
  gettimeofday(&begin_computepwgts, NULL);
  compute_pwgts<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_vwgt,graph->cuda_where,gainv);
  cudaDeviceSynchronize();
  gettimeofday(&end_computepwgts, NULL);
  computepwgts += (end_computepwgts.tv_sec - begin_computepwgts.tv_sec) * 1000 + (end_computepwgts.tv_usec - begin_computepwgts.tv_usec) / 1000.0;

  cudaDeviceSynchronize();
  gettimeofday(&begin_thrustreduce, NULL);
  for(int l = nvtxs;l != 1;l = (l + 512 - 1) / 512)
        reduction6<<<(l + 512 - 1) / 512,256>>>(gainv,l);
  cudaMemcpy(&graph->pwgts[0], &gainv[0],sizeof(int),cudaMemcpyDeviceToHost);
//   graph->pwgts[0] = thrust::reduce(thrust::device, gainv, gainv + nvtxs);
  graph->pwgts[1] = graph->tvwgt[0] - graph->pwgts[0];
  cudaDeviceSynchronize();
  gettimeofday(&end_thrustreduce, NULL);
  thrustreduce += (end_thrustreduce.tv_sec - begin_thrustreduce.tv_sec) * 1000 + (end_thrustreduce.tv_usec - begin_thrustreduce.tv_usec) / 1000.0;

  if((graph->pwgts[0] >= graph->tvwgt[0] * 0.5 / 1.03 && graph->pwgts[0] <= graph->tvwgt[0] * 0.5 * 1.03) && (graph->pwgts[1] >= graph->tvwgt[0] * 0.5 / 1.03 && graph->pwgts[1] <= graph->tvwgt[0] * 0.5 * 1.03)) /*printf("balance\n")*/;
  else
  {
    // printf("1\n");

    if(graph->pwgts[0] > graph->pwgts[1]) from = 0;
    else from = 1;
    to = (from + 1) % 2;

    cudaDeviceSynchronize();
    gettimeofday(&begin_computegain, NULL);
    compute_gainkp<<<(nvtxs + 127) / 128,128>>>(nvtxs, graph->nedges, from, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt, graph->cuda_where, kp);
    cudaDeviceSynchronize();
    gettimeofday(&end_computegain, NULL);
    computegain += (end_computegain.tv_sec - begin_computegain.tv_sec) * 1000 + (end_computegain.tv_usec - begin_computegain.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_thrustsort, NULL);
    thrust::sort(thrust::device, kp, kp + nvtxs, compRule());
    cudaDeviceSynchronize();
    gettimeofday(&end_thrustsort, NULL);
    thrustsort += (end_thrustsort.tv_sec - begin_thrustsort.tv_sec) * 1000 + (end_thrustsort.tv_usec - begin_thrustsort.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_computegainv, NULL);
    compute_gainvkp<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_vwgt,kp,gainv);
    cudaDeviceSynchronize();
    gettimeofday(&end_computegainv, NULL);
    computegainv += (end_computegainv.tv_sec - begin_computegainv.tv_sec) * 1000 + (end_computegainv.tv_usec - begin_computegainv.tv_usec) / 1000.0;

    cudaDeviceSynchronize();
    gettimeofday(&begin_inclusive, NULL);
    thrust::inclusive_scan(thrust::device, gainv, gainv + nvtxs, gainv);
    cudaDeviceSynchronize();
    gettimeofday(&end_inclusive, NULL);
    inclusive += (end_inclusive.tv_sec - begin_inclusive.tv_sec) * 1000 + (end_inclusive.tv_usec - begin_inclusive.tv_usec) / 1000.0;

    //move
    cudaDeviceSynchronize();
    gettimeofday(&begin_rebalance, NULL);
    rebalancekp<<<(nvtxs + 127) / 128,128>>>(nvtxs,(graph->pwgts[from] - graph->pwgts[to]) / 2,to,graph->cuda_where,kp,gainv);
    cudaDeviceSynchronize();
    gettimeofday(&end_rebalance, NULL);
    re_balance += (end_rebalance.tv_sec - begin_rebalance.tv_sec) * 1000 + (end_rebalance.tv_usec - begin_rebalance.tv_usec) / 1000.0;

  }
  
}

__global__ void select_bnd(int nvtxs, int *xadj, int *adjncy, int *where, int *bnd)
{
	int iii = threadIdx.x / 32;
	int ii  = blockIdx.x * 4 + iii;
    int tid = threadIdx.x % 32;

	// 启动�?1个Block�?128线程
	// 1个Warp处理一个顶点，将处理顶点的where存放于共享内�?
	// 可以减少32倍的where全局访问，但是会不会引起Bank冲突�?
	__shared__ int cache_twhere[4];		// 目的是减少global memory访问
	__shared__ int cache_flag[4];		// 目的是warp内通信

	if(tid == 0) cache_twhere[iii] = where[ii];
	if(tid == 1) cache_flag[iii]   = 0;
	__syncthreads();

	if(ii < nvtxs)
	{
		int i, j, wj, me, begin, end;

		me    = cache_twhere[iii];						// 会不会引起Bank冲突
		begin = xadj[ii];
		end   = xadj[ii + 1];

		for(i = begin + tid;i < end;i += 32)
		{
			j  = adjncy[i];
			wj = where[j];
			if(wj != me) cache_flag[iii]++;
			if(cache_flag[iii] != 0) break;
		}

		// 是边界顶点，bnd值为1；不是边界顶点，bnd值为1
		// 可以通过reduce计算边界顶点的数�?
		if(tid == 0)
		{
			if(cache_flag[iii] != 0) bnd[ii] = 1;
			else bnd[ii] = 0;
		}
	}
}

__global__ void compute_gain_shared(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *bnd, int *gain)
{
	int iii = threadIdx.x / 32;
	int ii  = blockIdx.x * 4 + iii;
	int tid = threadIdx.x % 32;

	// 启动�?1个Block�?128线程
	// 1个Warp处理一个顶点，将处理顶点的where存放于共享内�?
	// 可以减少32倍的where全局访问，但是会不会引起Bank冲突�?
	__shared__ int cache_twhere[4];		// 目的是减少global memory访问
	__shared__ int cache_ed[128];		// 目的是warp内reduce计算负责顶点的ed
	__shared__ int cache_id[128];		// 目的是warp内reduce计算负责顶点的id

	if(tid == 0) cache_twhere[iii] = where[ii];
	cache_ed[threadIdx.x] = 0;
	cache_id[threadIdx.x] = 0;
	__syncthreads();

	if(ii < nvtxs && bnd[ii] == 1)		//bnd[ii]理论上也可以像cache_twhere减少32倍的全局访问
	{
		int i, j, wj, me, begin ,end;
		me    = cache_twhere[iii];
		begin = xadj[ii];
		end   = xadj[ii + 1];

		for(i = begin;i < end;i += 32)
		{
			j  = adjncy[i];
			wj = where[j];
			if(me == wj) cache_id[threadIdx.x] += adjwgt[i];
			else cache_ed[threadIdx.x] += adjwgt[i];
		}

		// 对cache_id和cache_ed进行Warp内reduce
		if(tid < 16) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 16];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 16];
		}
		if(tid < 8) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 8];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 8];
		}
		if(tid < 4) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 4];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 4];
		}
		if(tid < 2) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 2];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 2];
		}
		if(tid < 1) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 1];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 1];
		}

		// 计算gain并存�?
		if(tid == 0) gain[ii] = cache_ed[threadIdx.x] - cache_id[threadIdx.x];
	}
}

__global__ void compute_infor_shared(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *gain)
{
	int iii = threadIdx.x / 32;
	int ii  = blockIdx.x * 4 + iii;
  	int tid = threadIdx.x % 32;

	// 启动�?1个Block�?128线程
	// 1个Warp处理一个顶点，将处理顶点的where存放于共享内�?
	// 可以减少32倍的where全局访问，但是会不会引起Bank冲突�?
	// __shared__ int cache_twhere[4];		// 目的是减少global memory访问
	__shared__ int cache_ed[128];		// 目的是warp内reduce计算负责顶点的ed
	__shared__ int cache_id[128];		// 目的是warp内reduce计算负责顶点的id

	// if(tid == 0) cache_twhere[iii] = where[ii];
	// cache_ed[threadIdx.x] = 0;
	// cache_id[threadIdx.x] = 0;
	// __syncthreads();

	if(ii < nvtxs)
	{
		// if(tid == 0) cache_twhere[iii] = where[ii];
		cache_ed[threadIdx.x] = 0;
		cache_id[threadIdx.x] = 0;
		__syncthreads();

		int i, j, wj, me, begin ,end, t;
		// me    = cache_twhere[iii];
		me    = where[ii];
		begin = xadj[ii];
		end   = xadj[ii + 1];

		for(i = begin;i < end;i += 32)
		{
			j  = adjncy[i];
			wj = where[j];
			if(me == wj) cache_id[threadIdx.x] += adjwgt[i];
			else cache_ed[threadIdx.x] += adjwgt[i];
		}

		// 对cache_id和cache_ed进行Warp内reduce
		if(tid < 16) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 16];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 16];
		}
		if(tid < 8) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 8];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 8];
		}
		if(tid < 4) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 4];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 4];
		}
		if(tid < 2) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 2];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 2];
		}
		if(tid < 1) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 1];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 1];
		}

		// 计算gain并存放（移动也可以放在这里）
		if(tid == 0) 
		{
			t = cache_ed[threadIdx.x] - cache_id[threadIdx.x];
			// printf("ii=%d gain=%d\n",ii,t);
			// t = 0;
			// gain[ii] = t;
			if(t > 0) where[ii] = (me + 1) % 2;
		}
	}
}

__global__ void compute_infor(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *gain)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, k, me, wk, begin, end, ed, id, t;
		me    = where[ii];
		begin = xadj[ii];
		end   = xadj[ii + 1];
		ed    = 0;
		id    = 0;

		for(i = begin;i < end;i++)
		{
			j  = adjncy[i];
			wk = where[j];

			if(me == wk) id += adjwgt[i];
			else ed += adjwgt[i];
		}
		
		t = ed - id;
		gain[ii] = t;
		if(t > 0) where[ii] = (me + 1) % 2;
	}
}

__global__ void compute_infor_atomic(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *gain)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, k, me, wk, begin, end, ed, id, t, p;
		p = 10;
		while(p >= 0)
		{
			me    = where[ii];
			begin = xadj[ii];
			end   = xadj[ii + 1];
			ed    = 0;
			id    = 0;

			for(i = begin;i < end;i++)
			{
				j  = adjncy[i];
				// wk = atomicAdd(&where[j], 0);
				wk = where[j];

				if(me == wk) id += adjwgt[i];
				else ed += adjwgt[i];
			}
			
			t = ed - id;
			// gain[ii] = t;
			if(t > 0) atomicExch(&where[ii], (me + 1) % 2);
			p --;
		}
	}
}

__global__ void greedy_refine(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *step)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, k, me, wk, begin, end, ed, id, t, p, s;
		
		s = nvtxs % 100;
		p = 0;

		while(p < s)
		{
			me    = where[ii];
			begin = xadj[ii];
			end   = xadj[ii + 1];
			ed    = 0;
			id    = 0;

			for(i = begin;i < end;i++)
			{
				j  = adjncy[i];
				wk = where[j];

				if(me == wk) id += adjwgt[i];
				else ed += adjwgt[i];
			}
			
			t = ed - id;

			if(t > 0)
			{
				if(atomicAdd(step, 1) <= p * 3)
				{
					atomicExch(&where[ii], (me + 1) % 2);
					break;
				}
			}
			
			p++;
		}
	}
}

__global__ void compute_balance(int nvtxs, int *vwgt, int *where, int *gain)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 0) gain[ii] = vwgt[ii];
		else gain[ii] = 0;
	}
}

__global__ void move_vertex(int nvtxs, int *where, int *gain)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(gain[ii] > 0) 
		{
			int me = where[ii];

			// 取余快还是if快？
			if(me == 0) where[ii] = 1;
			else where[ii] = 0;
		}
	}
}

__global__ void compute_infor_rebalance_shared(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, kp_t *gain_kp)
{
	int iii = threadIdx.x / 32;
	int ii  = blockIdx.x * 4 + iii;
  	int tid = threadIdx.x % 32;

	// 启动�?1个Block�?128线程
	// 1个Warp处理一个顶点，将处理顶点的where存放于共享内�?
	// 可以减少32倍的where全局访问，但是会不会引起Bank冲突�?
	// __shared__ int cache_twhere[4];		// 目的是减少global memory访问
	__shared__ int cache_ed[128];		// 目的是warp内reduce计算负责顶点的ed
	__shared__ int cache_id[128];		// 目的是warp内reduce计算负责顶点的id

	// if(tid == 0) cache_twhere[iii] = where[ii];
	// cache_ed[threadIdx.x] = 0;
	// cache_id[threadIdx.x] = 0;
	// __syncthreads();

	if(ii < nvtxs)
	{
		// if(tid == 0) cache_twhere[iii] = where[ii];
		cache_ed[threadIdx.x] = 0;
		cache_id[threadIdx.x] = 0;
		__syncthreads();

		int i, j, wj, me, begin ,end, t;
		// me    = cache_twhere[iii];
		me    = where[ii];
		begin = xadj[ii];
		end   = xadj[ii + 1];

		for(i = begin;i < end;i += 32)
		{
			j  = adjncy[i];
			wj = where[j];
			if(me == wj) cache_id[threadIdx.x] += adjwgt[i];
			else cache_ed[threadIdx.x] += adjwgt[i];
		}

		// 对cache_id和cache_ed进行Warp内reduce
		if(tid < 16) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 16];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 16];
		}
		if(tid < 8) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 8];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 8];
		}
		if(tid < 4) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 4];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 4];
		}
		if(tid < 2) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 2];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 2];
		}
		if(tid < 1) 
		{
			cache_ed[threadIdx.x] += cache_ed[threadIdx.x + 1];
			cache_id[threadIdx.x] += cache_id[threadIdx.x + 1];
		}

		// 计算gain并存放（移动也可以放在这里）
		if(tid == 0) 
		{
			gain_kp[ii].key = cache_ed[threadIdx.x] - cache_id[threadIdx.x];
			gain_kp[ii].ptr = ii;
		}
	}
}

__global__ void compute_infor_rebalance(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, kp_t *gain_kp)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int i, j, k, me, wk, begin, end, ed, id, t;
		me    = where[ii];
		begin = xadj[ii];
		end   = xadj[ii + 1];
		ed    = 0;
		id    = 0;

		for(i = begin;i < end;i++)
		{
			j  = adjncy[i];
			wk = where[j];

			if(me == wk) id += adjwgt[i];
			else ed += adjwgt[i];
		}

		gain_kp[ii].key = ed - id;
		gain_kp[ii].ptr = ii;
	}
}

__global__ void compute_vwgt(int nvtxs, int *vwgt, kp_t *gain_kp, int *gain)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int temp = gain_kp[ii].ptr;
		gain[ii] = vwgt[temp];
	}
}

__global__ void rebalance(int nvtxs, int to, int val, int *where, kp_t *gain_kp, int *gain)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int temp, ptr;
		temp = gain[ii];

		if(temp < val)
		{
			ptr = gain_kp[ii].key;
			where[ptr] = to;
		}
	}
}

//利用原子操作实现等待 向CPU串行方向靠拢
void FM_GPU(hunyuangraph_graph_t *graph)
{
	int nvtxs, nbnd, from, to, val;
	int *bnd, *gain, *moved, *step;
	kp_t *gain_kp;

	nvtxs = graph->nvtxs;

	cudaMalloc((void **)&gain,sizeof(int) * nvtxs);
	cudaMalloc((void **)&step,sizeof(int));

	// select_bnd和compute_gain_shared可以合并成一个计算全部顶点的ed，id，gain，bnd
	// 合并select_bnd和compute_gain_shared可做进一步优化，可以减少全局访存

	// 目前没有标记边界顶点（bnd�?
	// for(int i = 0;i < 5;i++)	// 将这个循环放入kernel
	// 	// compute_infor_shared<<<(nvtxs + 3) / 4, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
	// 		graph->cuda_where, gain);
		// compute_infor<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
		// 	graph->cuda_where, gain);
	// compute_infor_atomic<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
		graph->cuda_where, gain);

	greedy_refine<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
		graph->cuda_where, step);

	// 停止移动后观察是否平�?
	compute_balance<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_vwgt,graph->cuda_where,gain);

	graph->pwgts[0] = thrust::reduce(thrust::device, gain, gain + nvtxs);
	// for(int l = nvtxs;l != 1;l = (l + 512 - 1) / 512)
    //     reduction6<<<(l + 512 - 1) / 512,256>>>(gain,l);
	// cudaMemcpy(&graph->pwgts[0], &gain[0],sizeof(int),cudaMemcpyDeviceToHost);

	graph->pwgts[1] = graph->tvwgt[0] - graph->pwgts[0];

	// 若平衡，不操�?
	if((graph->pwgts[0] >= graph->tvwgt[0] * 0.5 / 1.03 && graph->pwgts[0] <= graph->tvwgt[0] * 0.5 * 1.03) && \
		(graph->pwgts[1] >= graph->tvwgt[0] * 0.5 / 1.03 && graph->pwgts[1] <= graph->tvwgt[0] * 0.5 * 1.03)) ;

	// 若不平衡则移动顶点至平衡
	else
	{
		printf("rebalance\n");

		// 从权重大的分区向权重小的分区移动顶点
		if(graph->pwgts[0] > graph->pwgts[1]) from = 0;
		else from = 1;
		to = (from + 1) % 2;

		// 计算当前分区状态的gain
		// compute_infor_rebalance_shared<<<(nvtxs + 3) / 4, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
			graph->cuda_where, gain_kp);
		cudaMalloc((void **)&gain_kp,sizeof(kp_t) * nvtxs);
		compute_infor_rebalance<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
			graph->cuda_where, gain_kp);

		// 按照gain的数值从大到小进行排�?
		thrust::sort(thrust::device, gain_kp, gain_kp + nvtxs, compRule());

		// 通过前缀和计算排序后顶点权重的和，据此判断哪些顶点移动就可以满足平衡
		compute_vwgt<<<(nvtxs + 127) / 128, 128>>>(nvtxs, graph->cuda_vwgt, gain_kp, gain);

		thrust::inclusive_scan(thrust::device, gain, gain + nvtxs, gain);

		// 移动顶点至平�?
		val = (graph->pwgts[from] - graph->pwgts[to]) / 2;
		rebalance<<<(nvtxs + 127) / 128, 128>>>(nvtxs, to, val, graph->cuda_where, gain_kp, gain);
	}

	cudaFree(gain);
	cudaFree(step);
	cudaFree(gain_kp);
}

__global__ void init_moved(int nvtxs, int *moved)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
		moved[ii] = -1;
}

__global__ void compute_gain(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *gain, int *moved)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs && moved[ii] == -1)
	{
		int i, j, k, me, wk, begin, end, ed, id, t;
		me    = where[ii];
		begin = xadj[ii];
		end   = xadj[ii + 1];
		ed    = 0;
		id    = 0;

		for(i = begin;i < end;i++)
		{
			j  = adjncy[i];
			wk = where[j];

			if(me == wk) id += adjwgt[i];
			else ed += adjwgt[i];
		}
		
		t = ed - id;
		gain[ii] = t;
	}
}

__global__ void judge_move(int id, int *mincut, int *newcut, int *max_element_ptr, int *tpwgts, int *pwgts, int *origdiff, \
	int *avgvwgt, int *mindiff, int *mincutorder, int *nswaps, int *limit, int *flag, int *where, int *moved, int *swaps)
{
	newcut[0] = mincut[0] - max_element_ptr[0];
	if((newcut[0] < mincut[0]&&abs(tpwgts[0]-pwgts[0])<=origdiff[0]+avgvwgt[0])|| 
        (newcut[0] == mincut[0]&&abs(tpwgts[0]-pwgts[0])<mindiff[0])){
        mincut[0]=newcut[0];
        mindiff[0]=abs(tpwgts[0]-pwgts[0]);
        mincutorder[0]=nswaps[0];

		// 
		int me = where[id];
		if(me == 0) where[id] = 1;
		else where[id] = 0;
		moved[id]=nswaps[0];
		swaps[nswaps[0]]=id;
    }
    else if(nswaps[0] - mincutorder[0] > limit[0]){ 
        flag[0] = 1;
    }
}

__global__ void move_max_index(int id, int *where, int *moved)
{
	int me = where[id];
	moved[id] = 1;
	if(me == 0) where[id] = 1;
	else where[id] = 0;
}

void Greedy_GPU(hunyuangraph_graph_t *graph, float *ntpwgts)
{
	int *gain, *moved, *swaps;
	int *max_element_ptr, max_index;

	int nvtxs  = graph->nvtxs;
	int *pwgts = graph->pwgts;

	int mincut, initcut, newcut;
	mincut = initcut = newcut = graph->mincut;

	int tpwgts[2];
	tpwgts[0]=graph->tvwgt[0]*ntpwgts[0];
  	tpwgts[1]=graph->tvwgt[0]-tpwgts[0];

	int limit=hunyuangraph_min(hunyuangraph_max(0.01*nvtxs,15),100);
  	int avgvwgt=hunyuangraph_min((pwgts[0]+pwgts[1])/20,2*(pwgts[0]+pwgts[1])/nvtxs);
	int origdiff=abs(tpwgts[0]-pwgts[0]);

	cudaMalloc((void **)&gain,sizeof(int) * nvtxs);
	cudaMalloc((void **)&moved,sizeof(int) * nvtxs);
	cudaMalloc((void **)&swaps,sizeof(int) * nvtxs);

	init_moved<<<(nvtxs + 127) / 128, 128>>>(nvtxs, moved);

	// 计算gain
	compute_gain<<<(nvtxs + 127) / 128, 128>>>(nvtxs, graph->cuda_xadj, graph->cuda_adjncy, graph->cuda_adjwgt,\
		graph->cuda_where, gain, moved);
	
	// 挑选最优顶�? 使用thrust::max_element获取最大值的迭代�?
	max_element_ptr = thrust::max_element(thrust::device, gain, gain + nvtxs);
	max_index =  max_element_ptr - gain;

	// 判断最优顶点是否移�?
	// judge_move<<<1,1>>>(max_index,);
	/* CPU判断规则
	if((newcut<mincut&&abs(tpwgts[0]-pwgts[0])<=origdiff+avgvwgt)|| 
        (newcut==mincut&&abs(tpwgts[0]-pwgts[0])<mindiff)){
        mincut=newcut;
        mindiff=abs(tpwgts[0]-pwgts[0]);
        mincutorder=nswaps;
    }
    else if(nswaps-mincutorder>limit){ 
        newcut+=(ed[higain]-id[higain]);
        hunyuangraph_add_sub(pwgts[from],pwgts[to],vwgt[higain]);
        break;
      }
	*/

	// 移动最优顶�?
	move_max_index<<<1,1>>>(max_index, graph->cuda_where, moved);
}

void FM_2WayCutRefine_GPU(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, float *ntpwgts)
{
	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_2way, NULL);
	// cudaMalloc((void**)&graph->cuda_xadj,sizeof(int) * (graph->nvtxs + 1));
	// cudaMalloc((void**)&graph->cuda_adjncy,sizeof(int) * graph->nedges);
	// cudaMalloc((void**)&graph->cuda_adjwgt,sizeof(int) * graph->nedges);
	// cudaMalloc((void**)&graph->cuda_vwgt,sizeof(int) * graph->nvtxs);
	cudaMalloc((void**)&graph->cuda_where,sizeof(int) * graph->nvtxs);

	// cudaMemcpy(graph->cuda_xadj,graph->xadj,sizeof(int) * (graph->nvtxs + 1), cudaMemcpyHostToDevice);
	// cudaMemcpy(graph->cuda_adjncy,graph->adjncy,sizeof(int) * graph->nedges, cudaMemcpyHostToDevice);
	// cudaMemcpy(graph->cuda_adjwgt,graph->adjwgt,sizeof(int) * graph->nedges, cudaMemcpyHostToDevice);
	// cudaMemcpy(graph->cuda_vwgt,graph->vwgt,sizeof(int) * graph->nvtxs, cudaMemcpyHostToDevice);
	cudaMemcpy(graph->cuda_where,graph->where,sizeof(int) * graph->nvtxs, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_2way, NULL);
	malloc_2way += (end_malloc_2way.tv_sec - begin_malloc_2way.tv_sec) * 1000 + (end_malloc_2way.tv_usec - begin_malloc_2way.tv_usec) / 1000.0;

	cudaDeviceSynchronize();
	gettimeofday(&begin_gpu_2way, NULL);
	// jetLP(graph);
	FM_GPU(graph);
	// Greedy_GPU(graph, ntpwgts);
	cudaDeviceSynchronize();
	gettimeofday(&end_gpu_2way, NULL);
	gpu_2way += (end_gpu_2way.tv_sec - begin_gpu_2way.tv_sec) * 1000 + (end_gpu_2way.tv_usec - begin_gpu_2way.tv_usec) / 1000.0;

	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_2way, NULL);
	cudaMemcpy(graph->where,graph->cuda_where,sizeof(int) * graph->nvtxs, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_2way, NULL);
	malloc_2way += (end_malloc_2way.tv_sec - begin_malloc_2way.tv_sec) * 1000 + (end_malloc_2way.tv_usec - begin_malloc_2way.tv_usec) / 1000.0;
}

/*Cpu graph refine two partitions*/
void hunyuangraph_cpu_2way_refine(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, float *ntpwgts, int iteration_num)
{
  int i,ii,j,k,kwgt,nvtxs,nbnd,nswaps,from,to,pass,limit,temp;
  int *xadj,*vwgt,*adjncy,*adjwgt,*where,*id,*ed,*bndptr,*bndlist,*pwgts;
  int *moved,*swaps,*perm;

  hunyuangraph_queue_t *queues[2];
  int higain,mincut, mindiff,origdiff,initcut,newcut,mincutorder,avgvwgt;
  int tpwgts[2];

  nvtxs=graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
  where=graph->where;
  id=graph->id;
  ed=graph->ed;
  pwgts=graph->pwgts;
  bndptr=graph->bndptr;
  bndlist=graph->bndlist;

  moved=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  swaps=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  perm=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  CPU_malloc(sizeof(int) * nvtxs * 3);

  tpwgts[0]=graph->tvwgt[0]*ntpwgts[0];
  tpwgts[1]=graph->tvwgt[0]-tpwgts[0];

  limit=hunyuangraph_min(hunyuangraph_max(0.01*nvtxs,15),100);
  avgvwgt=hunyuangraph_min((pwgts[0]+pwgts[1])/20,2*(pwgts[0]+pwgts[1])/nvtxs);

  queues[0]=hunyuangraph_queue_create(nvtxs);
  queues[1]=hunyuangraph_queue_create(nvtxs);

  origdiff=abs(tpwgts[0]-pwgts[0]);
  hunyuangraph_int_set_value(nvtxs,-1,moved);

  for(pass=0;pass<iteration_num;pass++){ 
    hunyuangraph_queue_reset(queues[0]);
    hunyuangraph_queue_reset(queues[1]);

    mincutorder=-1;
    newcut=mincut=initcut=graph->mincut;
    mindiff=abs(tpwgts[0]-pwgts[0]);
    nbnd=graph->nbnd;
    hunyuangraph_int_randarrayofp(nbnd,perm,nbnd,1); 

    for(ii=0;ii<nbnd;ii++){
      i=perm[ii];
      hunyuangraph_queue_insert(queues[where[bndlist[i]]],bndlist[i],ed[bndlist[i]]-id[bndlist[i]]);
    }       

    for(nswaps=0;nswaps<nvtxs;nswaps++){
      from=(tpwgts[0]-pwgts[0]<tpwgts[1]-pwgts[1]?0:1);
      to=(from+1)%2;

      if((higain=hunyuangraph_queue_top(queues[from]))==-1){
        break;
      }

      newcut-=(ed[higain]-id[higain]);
      hunyuangraph_add_sub(pwgts[to],pwgts[from],vwgt[higain]);

      if((newcut<mincut&&abs(tpwgts[0]-pwgts[0])<=origdiff+avgvwgt)|| 
          (newcut==mincut&&abs(tpwgts[0]-pwgts[0])<mindiff)){
        mincut=newcut;
        mindiff=abs(tpwgts[0]-pwgts[0]);
        mincutorder=nswaps;
      }
      else if(nswaps-mincutorder>limit){ 
        newcut+=(ed[higain]-id[higain]);
        hunyuangraph_add_sub(pwgts[from],pwgts[to],vwgt[higain]);
        break;
      }

      where[higain]=to;
      moved[higain]=nswaps;
      swaps[nswaps]=higain;

      hunyuangraph_swap(id[higain],ed[higain],temp);

      if(ed[higain]==0&&xadj[higain]<xadj[higain+1]){ 
        hunyuangraph_listdelete(nbnd,bndlist,bndptr,higain);
      }

      for(j=xadj[higain];j<xadj[higain+1];j++){
        k=adjncy[j];
        kwgt=(to==where[k]?adjwgt[j]:-adjwgt[j]);
        hunyuangraph_add_sub(id[k],ed[k],kwgt);

        if(bndptr[k]!=-1){ 
          if(ed[k]==0){ 
            hunyuangraph_listdelete(nbnd,bndlist,bndptr,k);
            
            if(moved[k]==-1){  
              hunyuangraph_queue_delete(queues[where[k]],k);
            }
          }
          else{ 
            if(moved[k]==-1){ 
              hunyuangraph_queue_update(queues[where[k]],k,ed[k]-id[k]);
            }
          }
        }
        else{
          if(ed[k]>0){  
            hunyuangraph_listinsert(nbnd,bndlist,bndptr,k);
            
            if(moved[k]==-1){ 
              hunyuangraph_queue_insert(queues[where[k]],k,ed[k]-id[k]);
            }
          }
        }
      }
    }

	// printf("2way_CPU:nvtxs=%d moved=%d\n",nvtxs,nswaps);

    for(i=0;i<nswaps;i++){
      moved[swaps[i]]=-1;  
    }

    for(nswaps--;nswaps>mincutorder;nswaps--){
      higain=swaps[nswaps];
      to=where[higain]=(where[higain]+1)%2;
      hunyuangraph_swap(id[higain],ed[higain],temp);

      if(ed[higain]==0&&bndptr[higain]!=-1&&xadj[higain]<xadj[higain+1]){
        hunyuangraph_listdelete(nbnd,bndlist,bndptr,higain);
      }
      else if(ed[higain]>0&&bndptr[higain]==-1){
        hunyuangraph_listinsert(nbnd,bndlist,bndptr,higain);
      }

      hunyuangraph_add_sub(pwgts[to],pwgts[(to+1)%2],vwgt[higain]);

      for(j=xadj[higain];j<xadj[higain+1];j++){
        k=adjncy[j];
        kwgt=(to==where[k]?adjwgt[j]:-adjwgt[j]);
        hunyuangraph_add_sub(id[k],ed[k],kwgt);

        if(bndptr[k]!=-1&&ed[k]==0){
          hunyuangraph_listdelete(nbnd,bndlist,bndptr,k);
        }
        if(bndptr[k]==-1&&ed[k]>0){
          hunyuangraph_listinsert(nbnd,bndlist,bndptr,k);
        }
      }
    }

    graph->mincut=mincut;
    graph->nbnd=nbnd;

    // printf("pass=%d nvtxs=%d\n",pass,nvtxs);
    // printf("graph->mincut=%d\n\n",graph->mincut);

    if(mincutorder<=0||mincut==initcut){
      break;
    }

  }

  hunyuangraph_queue_free(queues[0]);
  hunyuangraph_queue_free(queues[1]);

}

/*Cpu growbisection algorithm*/
void huyuangraph_cpu_growbisection(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, float *ntpwgts, int niparts)
{
  int i,j,k,nvtxs,dd,nleft,first,last,pwgts[2],oneminpwgt,onemaxpwgt, 
      bestcut=0,iter;

  int *xadj,*vwgt,*adjncy,*where;
  int *queue,*tra,*bestwhere;

  nvtxs=graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;

  hunyuangraph_allocate_cpu_2waymem(hunyuangraph_admin,graph);
  CPU_malloc(sizeof(int) * (nvtxs * 5 + 2));

  where=graph->where;

  bestwhere=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  queue=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  tra=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  CPU_malloc(sizeof(int) * nvtxs * 3);

  onemaxpwgt=hunyuangraph_admin->ubfactors[0]*graph->tvwgt[0]*ntpwgts[1];
  oneminpwgt=(1.0/hunyuangraph_admin->ubfactors[0])*graph->tvwgt[0]*ntpwgts[1]; 
  
  for (iter=0; iter<niparts; iter++){

    cudaDeviceSynchronize();
    gettimeofday(&begin_part_bfs,NULL);
    
    hunyuangraph_int_set_value(nvtxs,1,where);
    hunyuangraph_int_set_value(nvtxs,0,tra);

    pwgts[1]=graph->tvwgt[0];
    pwgts[0]=0;
    queue[0]=hunyuangraph_int_randinrange(nvtxs);
    tra[queue[0]]=1;
    first=0; 
    last=1;
    nleft=nvtxs-1;
    dd=0;

    for(;;){
      if(first==last){ 
        if(nleft==0||dd){
          break;
        }

        k=hunyuangraph_int_randinrange(nleft);

        for(i=0;i<nvtxs;i++){
          if(tra[i]==0){
            if(k==0){
              break;
            }
            else{
              k--;
            }
          }
        }

        queue[0]=i;
        tra[i]=1;
        first=0; 
        last=1;
        nleft--;
      }

      i=queue[first++];

      if(pwgts[0]>0&&pwgts[1]-vwgt[i]<oneminpwgt){
        dd=1;
        continue;
      }

      where[i]=0;

      hunyuangraph_add_sub(pwgts[0],pwgts[1],vwgt[i]);

      if(pwgts[1]<=onemaxpwgt){
        break;
      }

      dd=0;

      for(j=xadj[i];j<xadj[i+1];j++){
        k=adjncy[j];

        if(tra[k]==0){
          queue[last++]=k;
          tra[k]=1;
          nleft--;
        }
      }
    }

    cudaDeviceSynchronize();
    gettimeofday(&end_part_bfs,NULL);
    part_bfs += (end_part_bfs.tv_sec - begin_part_bfs.tv_sec) * 1000 + (end_part_bfs.tv_usec - begin_part_bfs.tv_usec) / 1000.0;

    hunyuangraph_compute_cpu_2wayparam(hunyuangraph_admin,graph);
    hunyuangraph_2way_bal(hunyuangraph_admin,graph,ntpwgts);

    cudaDeviceSynchronize();
    gettimeofday(&begin_part_2refine,NULL);
    hunyuangraph_cpu_2way_refine(hunyuangraph_admin,graph,ntpwgts,hunyuangraph_admin->iteration_num);
    cudaDeviceSynchronize();
    gettimeofday(&end_part_2refine,NULL);
    part_2refine += (end_part_2refine.tv_sec - begin_part_2refine.tv_sec) * 1000 + (end_part_2refine.tv_usec - begin_part_2refine.tv_usec) / 1000.0;
    
    if(iter==0||bestcut>graph->mincut){
      bestcut=graph->mincut;
      hunyuangraph_int_copy(nvtxs,where,bestwhere);
      
      if(bestcut==0){
        break;
      }
    }
  }

  graph->mincut=bestcut;
  hunyuangraph_int_copy(nvtxs,bestwhere,where);

}

/*Free graph params*/
void hunyuangraph_free_graph(hunyuangraph_graph_t **r_graph) 
{
  hunyuangraph_graph_t *graph;
  graph=*r_graph;

  free(graph->xadj);
  free(graph->vwgt);
  free(graph->adjncy);
  free(graph->adjwgt);
  free(graph->where);
  free(graph->pwgts);
  free(graph->id);
  free(graph->ed);
  free(graph->bndptr);
  free(graph->bndlist);
  free(graph->tvwgt);
  free(graph->tvwgt_reverse);
  free(graph->label);
  free(graph->cmap);
  free(graph);
  CPU_free(sizeof(graph->xadj) + sizeof(graph->vwgt) + sizeof(graph->adjncy) + sizeof(graph->adjwgt) + sizeof(graph->where) + \
  	sizeof(graph->pwgts) + sizeof(graph->id) + sizeof(graph->ed) + sizeof(graph->bndptr) + sizeof(graph->bndlist) + \
	sizeof(graph->tvwgt) + sizeof(graph->tvwgt_reverse) + sizeof(graph->label) + sizeof(graph->cmap) + sizeof(graph));

  *r_graph = NULL;
}

/*Cpu graph 2-way projection*/
void hunyuangraph_2way_project(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
  int i,j,istart,iend,nvtxs,nbnd,me,tid,ted;
  int *xadj,*adjncy,*adjwgt;
  int *cmap,*where,*bndptr,*bndlist;
  int *cwhere,*cbndptr;
  int *id,*ed;

  hunyuangraph_graph_t *cgraph;
  hunyuangraph_allocate_cpu_2waymem(hunyuangraph_admin,graph);
  CPU_malloc(sizeof(int) * (nvtxs * 5 + 2));

  cgraph=graph->coarser;
  cwhere=cgraph->where;
  cbndptr=cgraph->bndptr;
  nvtxs=graph->nvtxs;
  cmap=graph->cmap;
  xadj=graph->xadj;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
  where=graph->where;
  id=graph->id;
  ed=graph->ed;

  bndptr=hunyuangraph_int_set_value(nvtxs,-1,graph->bndptr);
  bndlist=graph->bndlist;

  for(i=0;i<nvtxs;i++){
    j=cmap[i];
    where[i]=cwhere[j];
    cmap[i]=cbndptr[j];
  }

	cudaDeviceSynchronize();
	gettimeofday(&begin_save_init, NULL);
  for(nbnd=0,i=0;i<nvtxs;i++){
    istart=xadj[i];
    iend=xadj[i+1];
    tid=ted=0;

    if(cmap[i]==-1){ 
      for(j=istart;j<iend;j++){
        tid+=adjwgt[j];
      }
    }
    else{ 
      me=where[i];

      for(j=istart;j<iend;j++){
        if(me==where[adjncy[j]]){
          tid += adjwgt[j];
        }
        else{
          ted+=adjwgt[j];
        }
      }
    }

    id[i]=tid;
    ed[i]=ted;

    if(ted>0||istart==iend){ 
      hunyuangraph_listinsert(nbnd,bndlist,bndptr,i);
    }

  }
  	cudaDeviceSynchronize();
	gettimeofday(&end_save_init, NULL);
	save_init += (end_save_init.tv_sec - begin_save_init.tv_sec) * 1000 + (end_save_init.tv_usec - begin_save_init.tv_usec) / 1000.0;

  graph->mincut=cgraph->mincut;
  graph->nbnd=nbnd;

  hunyuangraph_int_copy(2,cgraph->pwgts,graph->pwgts);
  hunyuangraph_free_graph(&graph->coarser);
  graph->coarser=NULL;

}

__global__ void projectback_init(int *where, int *cwhere, int *cmap, int nvtxs)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int t = cmap[ii];
		where[ii] = cwhere[t];
	}
}

int n = 0;

/*Cpu refinement algorithm*/
void hunyuangraph_cpu_refinement(hunyuangraph_admin_t *hunyuangraph_admin, \
hunyuangraph_graph_t *orggraph, hunyuangraph_graph_t *graph, float *tpwgts)
{
  hunyuangraph_compute_cpu_2wayparam(hunyuangraph_admin,graph);

  for(;;){
	hunyuangraph_2way_bal(hunyuangraph_admin,graph,tpwgts);
    // if(graph->nvtxs <= 10000) hunyuangraph_2way_bal(hunyuangraph_admin,graph,tpwgts);
    cudaDeviceSynchronize();
    gettimeofday(&begin_part_2refine,NULL);
	hunyuangraph_cpu_2way_refine(hunyuangraph_admin,graph,tpwgts,hunyuangraph_admin->iteration_num);
    /*if(graph->nvtxs <= 10000) hunyuangraph_cpu_2way_refine(hunyuangraph_admin,graph,tpwgts,hunyuangraph_admin->iteration_num); 
	else
	{
		FM_2WayCutRefine_GPU(hunyuangraph_admin, graph, tpwgts);
		n++;
	}*/
    cudaDeviceSynchronize();
    gettimeofday(&end_part_2refine,NULL);
    part_2refine += (end_part_2refine.tv_sec - begin_part_2refine.tv_sec) * 1000 + (end_part_2refine.tv_usec - begin_part_2refine.tv_usec) / 1000.0;
    
    if(graph==orggraph){
      	break;
    }

	// cudaFree(graph->cuda_xadj);
	// cudaFree(graph->cuda_adjncy);
	// cudaFree(graph->cuda_adjwgt);
	// cudaFree(graph->cuda_vwgt);
	// cudaFree(graph->cuda_where);

    graph=graph->finer;
	// printf("nvtxs=%d\n",graph->nvtxs);
    cudaDeviceSynchronize();
    gettimeofday(&begin_part_2map,NULL);
	hunyuangraph_2way_project(hunyuangraph_admin,graph);
    /*if(graph->nvtxs <= 10000) hunyuangraph_2way_project(hunyuangraph_admin,graph);
	else
	{
		int *cwhere;
		cudaMalloc((void**)&graph->cuda_where,sizeof(int) * graph->nvtxs);
		cudaMalloc((void**)&cwhere,sizeof(int) * graph->coarser->nvtxs);
		printf("a\n");
		cudaMemcpy(cwhere,graph->coarser->where,sizeof(int) * graph->nvtxs, cudaMemcpyDeviceToHost);
		printf("b\n");
		projectback_init<<<(graph->nvtxs + 127) / 128,128>>>(graph->cuda_where,cwhere,graph->cuda_cmap,graph->nvtxs);
		printf("c\n");
		graph->where = (int *)malloc(sizeof(int) * graph->nvtxs);
		graph->pwgts = (int *)malloc(sizeof(int) * 2);
		cudaMemcpy(graph->where,graph->cuda_where,sizeof(int) * graph->nvtxs, cudaMemcpyDeviceToHost);
		cudaFree(graph->cuda_where);
		cudaFree(graph->cuda_cmap);
		cudaFree(cwhere);
		printf("d\n");
	}*/
    cudaDeviceSynchronize();
    gettimeofday(&end_part_2map,NULL);
    part_2map += (end_part_2map.tv_sec - begin_part_2map.tv_sec) * 1000 + (end_part_2map.tv_usec - begin_part_2map.tv_usec) / 1000.0;
  }

}

/*************************************************************************/
/*! Computes the maximum load imbalance difference of a partitioning 
    solution over all the constraints. 
    The difference is defined with respect to the allowed maximum 
    unbalance for the respective constraint. 
 */
/**************************************************************************/ 
float ComputeLoadImbalanceDiff(hunyuangraph_graph_t *graph, int nparts, float *pijbm,float *ubvec)
{
	int  j, *pwgts;
	float max, cur;
	pwgts = graph->pwgts;
	max = -1.0;
	for (j=0; j<nparts; j++) 
	{
		cur = pwgts[j]*pijbm[j] - ubvec[0];
		if (cur > max)
			max = cur;
	}
	return max;
}

/*Cpu multilevel bisection algorithm*/
int hunyuangraph_cpu_mlevelbisect(hunyuangraph_admin_t *hunyuangraph_admin, \
hunyuangraph_graph_t *graph, float *tpwgts)
{
	int niparts,bestobj=0,curobj=0,*bestwhere=NULL;
	hunyuangraph_graph_t *tgraph = graph;
	hunyuangraph_graph_t *cgraph;
	double bestbal=0.0, curbal=0.0;

	hunyuangraph_compute_2way_balance(hunyuangraph_admin,graph,tpwgts);

	// printf("cnvtxs=%d cnedges=%d\n",graph->nvtxs,graph->nedges);
	/*if(graph->nvtxs > 10000) 
	{
		int nvtxs  = graph->nvtxs;
		int nedges = graph->nedges;

		cudaMalloc((void**)&graph->cuda_xadj,(nvtxs+1)*sizeof(int));
		cudaMalloc((void**)&graph->cuda_vwgt,nvtxs*sizeof(int));
		cudaMalloc((void**)&graph->cuda_adjncy,nedges*sizeof(int));
		cudaMalloc((void**)&graph->cuda_adjwgt,nedges*sizeof(int));

		cudaMemcpy(graph->cuda_xadj,graph->xadj,(nvtxs + 1) * sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(graph->cuda_adjncy,graph->adjncy,nedges*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(graph->cuda_adjwgt,graph->adjwgt,nedges*sizeof(int),cudaMemcpyHostToDevice);
		cudaMemcpy(graph->cuda_vwgt,graph->vwgt,nvtxs*sizeof(int),cudaMemcpyHostToDevice);

		int level = 0;
		hunyuangraph_admin->maxvwgt = 1.5 * graph->tvwgt[0] / hunyuangraph_admin->Coarsen_threshold; 
		do
		{
			hunyuangraph_malloc_coarseninfo(hunyuangraph_admin,graph);
		
			hunyuangraph_gpu_match(hunyuangraph_admin,graph);

			graph->cmap=(int*)malloc(sizeof(int)*(graph->nvtxs));
			cudaMemcpy(graph->cmap, graph->cuda_cmap, graph->nvtxs * sizeof(int), cudaMemcpyDeviceToHost);

			// cudaFree(graph->cuda_xadj);
			// cudaFree(graph->cuda_vwgt);
			// cudaFree(graph->cuda_adjncy);
			// cudaFree(graph->cuda_adjwgt);
			// cudaFree(graph->cuda_cmap);

			graph = graph->coarser;

			cudaDeviceSynchronize();
			gettimeofday(&begin_save_init, NULL);
			graph->xadj   = (int *)malloc(sizeof(int) * (graph->nvtxs + 1)); 
			graph->vwgt   = (int *)malloc(sizeof(int) * graph->nvtxs); 
			graph->adjncy = (int *)malloc(sizeof(int) * graph->nedges);
			graph->adjwgt = (int *)malloc(sizeof(int) * graph->nedges);
			cudaMemcpy(graph->xadj, graph->cuda_xadj, (graph->nvtxs + 1) * sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(graph->vwgt, graph->cuda_vwgt, graph->nvtxs * sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(graph->adjncy, graph->cuda_adjncy, graph->nedges * sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(graph->adjwgt, graph->cuda_adjwgt, graph->nedges * sizeof(int), cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			gettimeofday(&end_save_init, NULL);
			save_init += (end_save_init.tv_sec - begin_save_init.tv_sec) * 1000 + (end_save_init.tv_usec - begin_save_init.tv_usec) / 1000.0;

			level++;
			// printf("level=%d\n",level);
		}while(
			graph->nvtxs > 10000 && \
			graph->nvtxs > hunyuangraph_admin->Coarsen_threshold && \
			graph->nvtxs < 0.85 * graph->finer->nvtxs && \
			graph->nedges > graph->nvtxs / 2); 

		// graph->xadj   = (int *)malloc(sizeof(int) * (graph->nvtxs + 1)); 
		// graph->vwgt   = (int *)malloc(sizeof(int) * graph->nvtxs); 
		// graph->adjncy = (int *)malloc(sizeof(int) * graph->nedges);
		// graph->adjwgt = (int *)malloc(sizeof(int) * graph->nedges);
		// cudaMemcpy(graph->xadj, graph->cuda_xadj, (graph->nvtxs + 1) * sizeof(int), cudaMemcpyDeviceToHost);
		// cudaMemcpy(graph->vwgt, graph->cuda_vwgt, graph->nvtxs * sizeof(int), cudaMemcpyDeviceToHost);
		// cudaMemcpy(graph->adjncy, graph->cuda_adjncy, graph->nedges * sizeof(int), cudaMemcpyDeviceToHost);
		// cudaMemcpy(graph->adjwgt, graph->cuda_adjwgt, graph->nedges * sizeof(int), cudaMemcpyDeviceToHost);

		cudaFree(graph->cuda_xadj);
		cudaFree(graph->cuda_vwgt);
		cudaFree(graph->cuda_adjncy);
		cudaFree(graph->cuda_adjwgt);
	}*/
	/*cgraph=hunyuangraph_cpu_coarsen(hunyuangraph_admin,graph);
	// printf("cgraph->nvtxs=%d\n",cgraph->nvtxs);
	niparts=5;

	huyuangraph_cpu_growbisection(hunyuangraph_admin,cgraph,tpwgts,niparts);

	graph = tgraph;
	hunyuangraph_cpu_refinement(hunyuangraph_admin,graph,cgraph,tpwgts);

	curobj=graph->mincut;
	bestobj=curobj;

	if(bestobj!=curobj){
		hunyuangraph_int_copy(graph->nvtxs,bestwhere,graph->where);
		hunyuangraph_compute_cpu_2wayparam(hunyuangraph_admin,graph);
	}*/

	if (hunyuangraph_admin->ncuts > 1)
	{
    	bestwhere = (int *)malloc(sizeof(int) * graph->nvtxs);
		CPU_malloc(sizeof(int) * graph->nvtxs);
	}

  	for (int i = 0; i < hunyuangraph_admin->ncuts; i++) 
	{
		cgraph=hunyuangraph_cpu_coarsen(hunyuangraph_admin,graph);

		niparts = (cgraph->nvtxs <= hunyuangraph_admin->Coarsen_threshold ? 5 : 7);
		huyuangraph_cpu_growbisection(hunyuangraph_admin,cgraph,tpwgts,niparts);

		hunyuangraph_cpu_refinement(hunyuangraph_admin,graph,cgraph,tpwgts);

		curobj = graph->mincut;
		curbal = ComputeLoadImbalanceDiff(graph, 2, hunyuangraph_admin->part_balance, hunyuangraph_admin->ubfactors);

		if (i == 0  || (curbal <= 0.0005 && bestobj > curobj) || (bestbal > 0.0005 && curbal < bestbal)) 
		{
			bestobj = curobj;
			bestbal = curbal;
			if (i < hunyuangraph_admin->ncuts-1)
				hunyuangraph_int_copy(graph->nvtxs, graph->where, bestwhere);
    	}

		if (bestobj == 0)
			break;

		if (i < hunyuangraph_admin->ncuts-1)
			hunyuangraph_free_graph(&graph);
  	}

	if (bestobj != curobj) {
		hunyuangraph_int_copy(graph->nvtxs, bestwhere, graph->where);
		hunyuangraph_compute_cpu_2wayparam(hunyuangraph_admin,graph);
	}

  	return bestobj;
}

/*Set split graph params*/
hunyuangraph_graph_t *hunyuangraph_set_splitgraph(hunyuangraph_graph_t *graph, int snvtxs, int snedges)
{
  hunyuangraph_graph_t *sgraph;
  sgraph=hunyuangraph_create_cpu_graph();

  sgraph->nvtxs=snvtxs;
  sgraph->nedges=snedges;

  sgraph->xadj=(int*)malloc(sizeof(int)*(snvtxs+1));
  sgraph->vwgt=(int*)malloc(sizeof(int)*(snvtxs+1));
  sgraph->adjncy=(int*)malloc(sizeof(int)*(snedges));
  sgraph->adjwgt=(int*)malloc(sizeof(int)*(snedges));
  sgraph->label=(int*)malloc(sizeof(int)*(snvtxs));
  sgraph->tvwgt=(int*)malloc(sizeof(int));
  sgraph->tvwgt_reverse=(float*)malloc(sizeof(float));

  return sgraph;

}

/*Split graph to lgraph and rgraph*/
void hunyuangraph_splitgraph(hunyuangraph_admin_t *hunyuangraph_admin, \
hunyuangraph_graph_t *graph, hunyuangraph_graph_t **r_lgraph, hunyuangraph_graph_t **r_rgraph)
{
  int i,j,k,l,istart,iend,mypart,nvtxs,snvtxs[2],snedges[2];
  int *xadj,*vwgt,*adjncy,*adjwgt,*label,*where,*bndptr;
  int *sxadj[2],*svwgt[2],*sadjncy[2],*sadjwgt[2],*slabel[2];
  int *rename;
  int *temp_adjncy,*temp_adjwgt;

  hunyuangraph_graph_t *lgraph,*rgraph;

  nvtxs=graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
  label=graph->label;
  where=graph->where;
  bndptr=graph->bndptr;

  rename=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  CPU_malloc(sizeof(int) * nvtxs);
  snvtxs[0]=snvtxs[1]=snedges[0]=snedges[1]=0;

  for(i=0;i<nvtxs;i++){
    k=where[i];
    rename[i]=snvtxs[k]++;
    snedges[k]+=xadj[i+1]-xadj[i];
  }

  lgraph=hunyuangraph_set_splitgraph(graph,snvtxs[0],snedges[0]);
  CPU_malloc(sizeof(int) * (snedges[0] * 2 + snvtxs[0] * 3 + 3));
  sxadj[0]=lgraph->xadj;
  svwgt[0]=lgraph->vwgt;
  sadjncy[0]=lgraph->adjncy; 	
  sadjwgt[0]=lgraph->adjwgt; 
  slabel[0]=lgraph->label;

  rgraph=hunyuangraph_set_splitgraph(graph,snvtxs[1],snedges[1]);
  CPU_malloc(sizeof(int) * (snedges[1] * 2 + snvtxs[1] * 3 + 3));
  sxadj[1]=rgraph->xadj;
  svwgt[1]=rgraph->vwgt;
  sadjncy[1]=rgraph->adjncy; 	
  sadjwgt[1]=rgraph->adjwgt; 
  slabel[1]=rgraph->label;

  snvtxs[0]=snvtxs[1]=snedges[0]=snedges[1]=0;
  sxadj[0][0]=sxadj[1][0]=0;

  for(i=0;i<nvtxs;i++){
    mypart=where[i];
    istart=xadj[i];
    iend=xadj[i+1];

    if(bndptr[i]==-1){ 
      temp_adjncy=sadjncy[mypart]+snedges[mypart]-istart;
      temp_adjwgt=sadjwgt[mypart]+snedges[mypart]-istart;

      for(j=istart;j<iend;j++){
        temp_adjncy[j]=adjncy[j];
        temp_adjwgt[j]=adjwgt[j]; 
      }

      snedges[mypart]+=iend-istart;
    }
    else{
      temp_adjncy=sadjncy[mypart];
      temp_adjwgt=sadjwgt[mypart];
      l=snedges[mypart];

      for(j=istart;j<iend;j++){
        k=adjncy[j];
        
        if(where[k]==mypart){
          temp_adjncy[l]=k;
          temp_adjwgt[l++]=adjwgt[j]; 
        }
      }
      snedges[mypart]=l;
    }

    svwgt[mypart][snvtxs[mypart]]=vwgt[i];
    // printf("i=%d label[i]=%d\n",i,label[i]);
    slabel[mypart][snvtxs[mypart]]=label[i];
    sxadj[mypart][++snvtxs[mypart]]=snedges[mypart];
  }

  for(mypart=0;mypart<2;mypart++){
    iend=sxadj[mypart][snvtxs[mypart]];
    temp_adjncy=sadjncy[mypart];

    for(i=0;i<iend;i++){ 
      temp_adjncy[i]=rename[temp_adjncy[i]];
    }
  }

  lgraph->nedges=snedges[0];
  rgraph->nedges=snedges[1];

//   printf("CPU:lnedges=%d rnedges=%d nedges=%d +:%d\n",lgraph->nedges,rgraph->nedges,graph->nedges,lgraph->nedges + rgraph->nedges + graph->mincut);
	/*printf("CPU_lgraph:\n");
  	printf("lgraph->nvtxs:%d lgraph->nedges:%d\n",lgraph->nvtxs,lgraph->nedges);
	printf("lgraph->vwgt:\n");
  	for(i = 0;i < lgraph->nvtxs;i++)
	{
		printf("%d ",lgraph->vwgt[i]);
	}
	printf("\nlgraph->label:\n");
  	for(i = 0;i < lgraph->nvtxs;i++)
	{
		printf("%d ",lgraph->label[i]);
	}
	printf("\nlgraph->xadj:\n");
	for(i = 0;i <= lgraph->nvtxs;i++)
	{
		printf("%d ",lgraph->xadj[i]);
	}
	printf("\nlgraph->adjncy:\n");
	for(i = 0;i < lgraph->nvtxs;i++)
	{
		for(int j = lgraph->xadj[i];j < lgraph->xadj[i + 1];j++)
		{
			printf("%d ",lgraph->adjncy[j]);
		}
		printf("\n");
	}
	printf("CPU_rgraph:\n");
	printf("rgraph->nvtxs:%d rgraph->nedges:%d\n",rgraph->nvtxs,rgraph->nedges);
	printf("rgraph->vwgt:\n");
  	for(i = 0;i < rgraph->nvtxs;i++)
	{
		printf("%d ",rgraph->vwgt[i]);
	}
	printf("\ngraph->label:\n");
  	for(i = 0;i < rgraph->nvtxs;i++)
	{
		printf("%d ",rgraph->label[i]);
	}
	printf("\nrgraph->xadj:\n");
	for(i = 0;i <= rgraph->nvtxs;i++)
	{
		printf("%d ",rgraph->xadj[i]);
	}
	printf("\nrgraph->adjncy:\n");
	for(i = 0;i < rgraph->nvtxs;i++)
	{
		for(int j = rgraph->xadj[i];j < rgraph->xadj[i + 1];j++)
		{
			printf("%d ",rgraph->adjncy[j]);
		}
		printf("\n");
	}*/

  hunyuangraph_set_graph_tvwgt(lgraph);
  hunyuangraph_set_graph_tvwgt(rgraph);

	// printf("lgraph->tvwgt=%d lgraph->tvwgt_reverse=%f rgraph->tvwgt=%d rgraph->tvwgt_reverse=%f\n",lgraph->tvwgt[0],lgraph->tvwgt_reverse[0],rgraph->tvwgt[0],rgraph->tvwgt_reverse[0]);


  *r_lgraph=lgraph;
  *r_rgraph=rgraph;

}

void hunyuangraph_splitgraph_first(hunyuangraph_admin_t *hunyuangraph_admin, \
	hunyuangraph_graph_t *graph, hunyuangraph_graph_t **r_lgraph, hunyuangraph_graph_t **r_rgraph)
{
  int i,j,k,l,istart,iend,mypart,nvtxs,snvtxs[2],snedges[2];
  int *xadj,*vwgt,*adjncy,*adjwgt,*label,*where,*bndptr;
  int *sxadj[2],*svwgt[2],*sadjncy[2],*sadjwgt[2],*slabel[2];
  int *rename;
  int *temp_adjncy,*temp_adjwgt;

  hunyuangraph_graph_t *lgraph,*rgraph;

  nvtxs=graph->nvtxs;
  xadj=graph->xadj;
  vwgt=graph->vwgt;
  adjncy=graph->adjncy;
  adjwgt=graph->adjwgt;
//   label=graph->label;
  where=graph->where;
  bndptr=graph->bndptr;

  rename=hunyuangraph_int_malloc_space(hunyuangraph_admin,nvtxs);
  CPU_malloc(sizeof(int) * nvtxs);
  snvtxs[0]=snvtxs[1]=snedges[0]=snedges[1]=0;

  for(i=0;i<nvtxs;i++){
    k=where[i];
    rename[i]=snvtxs[k]++;
    snedges[k]+=xadj[i+1]-xadj[i];
  }

  lgraph=hunyuangraph_set_splitgraph(graph,snvtxs[0],snedges[0]);
  CPU_malloc(sizeof(int) * (snedges[0] * 2 + snvtxs[0] * 3 + 3));
  sxadj[0]=lgraph->xadj;
  svwgt[0]=lgraph->vwgt;
  sadjncy[0]=lgraph->adjncy; 	
  sadjwgt[0]=lgraph->adjwgt; 
  slabel[0]=lgraph->label;

  rgraph=hunyuangraph_set_splitgraph(graph,snvtxs[1],snedges[1]);
  CPU_malloc(sizeof(int) * (snedges[1] * 2 + snvtxs[1] * 3 + 3));
  sxadj[1]=rgraph->xadj;
  svwgt[1]=rgraph->vwgt;
  sadjncy[1]=rgraph->adjncy; 	
  sadjwgt[1]=rgraph->adjwgt; 
  slabel[1]=rgraph->label;

  snvtxs[0]=snvtxs[1]=snedges[0]=snedges[1]=0;
  sxadj[0][0]=sxadj[1][0]=0;

  for(i=0;i<nvtxs;i++){
    mypart=where[i];
    istart=xadj[i];
    iend=xadj[i+1];

    if(bndptr[i]==-1){ 
      temp_adjncy=sadjncy[mypart]+snedges[mypart]-istart;
      temp_adjwgt=sadjwgt[mypart]+snedges[mypart]-istart;

      for(j=istart;j<iend;j++){
        temp_adjncy[j]=adjncy[j];
        temp_adjwgt[j]=adjwgt[j]; 
      }

      snedges[mypart]+=iend-istart;
    }
    else{
      temp_adjncy=sadjncy[mypart];
      temp_adjwgt=sadjwgt[mypart];
      l=snedges[mypart];

      for(j=istart;j<iend;j++){
        k=adjncy[j];
        
        if(where[k]==mypart){
          temp_adjncy[l]=k;
          temp_adjwgt[l++]=adjwgt[j]; 
        }
      }
      snedges[mypart]=l;
    }

    svwgt[mypart][snvtxs[mypart]]=vwgt[i];
    slabel[mypart][snvtxs[mypart]]=i;
    sxadj[mypart][++snvtxs[mypart]]=snedges[mypart];
  }

  for(mypart=0;mypart<2;mypart++){
    iend=sxadj[mypart][snvtxs[mypart]];
    temp_adjncy=sadjncy[mypart];

    for(i=0;i<iend;i++){ 
      temp_adjncy[i]=rename[temp_adjncy[i]];
    }
  }

  lgraph->nedges=snedges[0];
  rgraph->nedges=snedges[1];

  	/*printf("CPU_graph:\n");
  	printf("graph->nvtxs:%d graph->nedges:%d\n",graph->nvtxs,graph->nedges);
	printf("\ngraph->xadj:\n");
	for(i = 0;i < graph->nvtxs;i++)
	{
		printf("i=%d where[i]=%d rename=%d\n",i,graph->where[i],rename[i]);
	}
	for(i = 0;i <= graph->nvtxs;i++)
	{
		printf("%d ",graph->xadj[i]);
	}
	printf("\ngraph->adjncy:\n");
	for(i = 0;i < graph->nvtxs;i++)
	{
		for(int j = graph->xadj[i];j < graph->xadj[i + 1];j++)
		{
			printf("i=%d j=%d where[i]=%d where[j]=%d\n",i,graph->adjncy[j],graph->where[i],graph->where[graph->adjncy[j]]);
		}
	}*/

	/*printf("CPU_lgraph:\n");
  	printf("lgraph->nvtxs:%d lgraph->nedges:%d\n",lgraph->nvtxs,lgraph->nedges);
	printf("lgraph->vwgt:\n");
  	for(i = 0;i < lgraph->nvtxs;i++)
	{
		printf("%d ",lgraph->vwgt[i]);
	}
	printf("\nlgraph->label:\n");
  	for(i = 0;i < lgraph->nvtxs;i++)
	{
		printf("%d ",lgraph->label[i]);
	}
	printf("\nlgraph->xadj:\n");
	for(i = 0;i <= lgraph->nvtxs;i++)
	{
		printf("%d ",lgraph->xadj[i]);
	}
	printf("\nlgraph->adjncy:\n");
	for(i = 0;i < lgraph->nvtxs;i++)
	{
		for(int j = lgraph->xadj[i];j < lgraph->xadj[i + 1];j++)
		{
			printf("%d ",lgraph->adjncy[j]);
		}
		printf("\n");
	}
	printf("CPU_rgraph:\n");
	printf("rgraph->nvtxs:%d rgraph->nedges:%d\n",rgraph->nvtxs,rgraph->nedges);
	printf("rgraph->vwgt:\n");
  	for(i = 0;i < rgraph->nvtxs;i++)
	{
		printf("%d ",rgraph->vwgt[i]);
	}
	printf("\ngraph->label:\n");
  	for(i = 0;i < rgraph->nvtxs;i++)
	{
		printf("%d ",rgraph->label[i]);
	}
	printf("\nrgraph->xadj:\n");
	for(i = 0;i <= rgraph->nvtxs;i++)
	{
		printf("%d ",rgraph->xadj[i]);
	}
	printf("\nrgraph->adjncy:\n");
	for(i = 0;i < rgraph->nvtxs;i++)
	{
		for(int j = rgraph->xadj[i];j < rgraph->xadj[i + 1];j++)
		{
			printf("%d ",rgraph->adjncy[j]);
		}
		printf("\n");
	}*/

  hunyuangraph_set_graph_tvwgt(lgraph);
  hunyuangraph_set_graph_tvwgt(rgraph);

//   printf("lgraph->tvwgt=%d lgraph->tvwgt_reverse=%f rgraph->tvwgt=%d rgraph->tvwgt_reverse=%f\n",lgraph->tvwgt[0],lgraph->tvwgt_reverse[0],rgraph->tvwgt[0],rgraph->tvwgt_reverse[0]);


  *r_lgraph=lgraph;
  *r_rgraph=rgraph;

}

__global__ void compute_lnedges(int nvtxs, int *xadj, int *adjncy, int *where, int *ltxadj)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 0)
		{
			int i, j, l, begin, end;
			begin = xadj[ii];
			end   = xadj[ii + 1];
			l     = 0;

			for(i = begin;i < end;i++)
			{
				j = adjncy[i];
				if(where[j] == 0) l ++;
			}
			ltxadj[ii] = l;
		}
		else
		{
			ltxadj[ii] = 0;
		}
	}
}

__global__ void compute_rnedges(int nvtxs, int *xadj, int *adjncy, int *where, int *rtxadj)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		// printf("ii=%d 1\n",ii);
		if(where[ii] == 1)
		{
			// printf("ii=%d 2\n",ii);
			int i, j, l, begin, end;
			begin = xadj[ii];
			end   = xadj[ii + 1];
			l     = 0;

			// printf("ii=%d 3\n",ii);
			for(i = begin;i < end;i++)
			{
				j = adjncy[i];
				if(where[j] == 1) l ++;
			}
			rtxadj[ii] = l;
		}
		else
		{
			rtxadj[ii] = 0;
		}
		// printf("ii=%d rtxadj=%d\n",ii,rtxadj[ii]);
	}
}

__global__ void compute_lmap(int nvtxs, int *where, int *lmap)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 0) lmap[ii] = 1;
		else lmap[ii] = 0;
	}
}

__global__ void compute_rmap(int nvtxs, int *where, int *rmap)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 1) rmap[ii] = 1;
		else rmap[ii] = 0;
	}
}

__global__ void set_lxadj(int nvtxs, int lnvtxs, int lnedges, int *xadj, int *where, int *lmap, int *ltxadj, int *lxadj)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 0) 
		{
			int val, ptr;
			val = ltxadj[ii];
			ptr = lmap[ii] - 1;
			lxadj[ptr] = val;
		}
	}
	else if(ii == nvtxs) lxadj[lnvtxs] = lnedges;
}

__global__ void set_rxadj(int nvtxs, int rnvtxs, int rnedges, int *xadj, int *where, int *rmap, int *rtxadj, int *rxadj)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 1) 
		{
			int val, ptr;
			val = rtxadj[ii];
			ptr = rmap[ii] - 1;
			rxadj[ptr] = val;
		}
	}
	else if(ii == nvtxs) rxadj[rnvtxs] = rnedges;
}

__global__ void set_ladjncy_ladjwgt(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *lxadj, int *ladjncy,\
	int *ladjwgt, int *lmap)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 0)
		{
			int i, j, begin, end, ptr;
			begin = xadj[ii];
			end   = xadj[ii + 1];
			ptr   = lmap[ii] - 1;
			ptr   = lxadj[ptr];

			for(i = begin;i < end;i++)
			{
				j = adjncy[i];
				if(where[j] == 0) 
				{
					ladjncy[ptr] = lmap[j] - 1;
					ladjwgt[ptr] = adjwgt[i];
					ptr++;
				}
			}
		}
	}
}

__global__ void set_radjncy_radjwgt(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where, int *rxadj, int *radjncy,\
	int *radjwgt, int *rmap)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(where[ii] == 1)
		{
			int i, j, begin, end, ptr;
			begin = xadj[ii];
			end   = xadj[ii + 1];
			ptr   = rmap[ii] - 1;
			ptr   = rxadj[ptr];

			for(i = begin;i < end;i++)
			{
				j = adjncy[i];
				if(where[j] == 1) 
				{
					radjncy[ptr] = rmap[j] - 1;
					radjwgt[ptr] = adjwgt[i];
					ptr++;
				}
			}
		}
	}
}

__global__ void set_lvwgt(int nvtxs, int *where, int *vwgt, int *lmap, int *lvwgt)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs && where[ii] == 0)
	{
		int ptr = lmap[ii] - 1;
		lvwgt[ptr] = vwgt[ii];
	}
}

__global__ void set_rvwgt(int nvtxs, int *where, int *vwgt, int *rmap, int *rvwgt)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs && where[ii] == 1)
	{
		int ptr = rmap[ii] - 1;
		rvwgt[ptr] = vwgt[ii];
	}
}

__global__ void set_llabel0(int nvtxs, int *where, int *lmap, int *llabel)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs && where[ii] == 0)
	{
		int ptr = lmap[ii] - 1;
		llabel[ptr] = ii;
	}
}

__global__ void set_rlabel0(int nvtxs, int *where, int *lmap, int *rlabel)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs && where[ii] == 1)
	{
		int ptr = lmap[ii] - 1;
		rlabel[ptr] = ii;
	}
}

__global__ void set_llabel1(int nvtxs, int *where, int *lmap, int *label, int *llabel)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs && where[ii] == 0)
	{
		int ptr = lmap[ii] - 1;
		llabel[ptr] = label[ii];
	}
}

__global__ void set_rlabel1(int nvtxs, int *where, int *lmap, int *label, int *rlabel)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs && where[ii] == 1)
	{
		int ptr = lmap[ii] - 1;
		rlabel[ptr] = label[ii];
	}
}

void set_subgraph_tvwgt(hunyuangraph_graph_t *graph, hunyuangraph_graph_t *lgraph, hunyuangraph_graph_t *rgraph)
{
	if(lgraph->tvwgt == NULL)
        lgraph->tvwgt = (int*)malloc(sizeof(int));
    if(lgraph->tvwgt_reverse == NULL)
        lgraph->tvwgt_reverse = (float*)malloc(sizeof(float));

	if(rgraph->tvwgt == NULL)
        rgraph->tvwgt = (int*)malloc(sizeof(int));
    if(rgraph->tvwgt_reverse == NULL)
        rgraph->tvwgt_reverse = (float*)malloc(sizeof(float));

	lgraph->tvwgt[0] = graph->pwgts[0];
	lgraph->tvwgt_reverse[0] = 1.0 / (lgraph->tvwgt[0] > 0 ? lgraph->tvwgt[0] : 1);

	rgraph->tvwgt[0] = graph->pwgts[1];
	rgraph->tvwgt_reverse[0] = 1.0 / (rgraph->tvwgt[0] > 0 ? rgraph->tvwgt[0] : 1);
}

__global__ void exam_temp_scan(int nvtxs, int *xadj, int *adjncy, int *where, int *temp_scan)
{
	for(int i = 0;i < nvtxs;i++)
	{
		for(int j = xadj[i];j < xadj[i + 1];j++)
			printf("%d ",adjncy[j]);
		printf("\n");
		for(int j = xadj[i];j < xadj[i + 1];j++)
			printf("%d ",where[adjncy[j]]);
		printf("\n");
		for(int j = xadj[i];j < xadj[i + 1];j++)
			printf("%d ",temp_scan[j]);
		printf("\n");
	}
}

__global__ void exam_lmap(int nvtxs, int *where, int *lmap)
{
	for(int i = 0;i < nvtxs;i++)
	{
		printf("ii=%d where=%d lmap=%d\n",i,where[i],lmap[i] - 1);
	}
}

__global__ void set_where(int *where)
{
	where[0] = 0;
	where[1] = 0;
	where[2] = 1;
	where[3] = 0;
	where[4] = 1;
}

__global__ void exam_graph(int nvtxs, int nedges, int *xadj, int *adjncy, int *vwgt, int *label)
{
	printf("lgraph->nvtxs:%d lgraph->nedges:%d\n",nvtxs,nedges);
  	printf("lgraph->vwgt:\n");
	for(int i = 0;i < nvtxs;i++)
	{
		printf("%d ",vwgt[i]);
	}
	printf("\nlgraph->label:\n");
	for(int i = 0;i < nvtxs;i++)
	{
		printf("%d ",label[i]);
	}
	printf("\nlgraph->xadj:\n");
	for(int i = 0;i <= nvtxs;i++)
	{
		printf("%d ",xadj[i]);
	}
	printf("\nlgraph->adjncy:\n");
	for(int i = 0;i < nvtxs;i++)
	{
		for(int j = xadj[i];j < xadj[i + 1];j++)
		{
			printf("%d ",adjncy[j]);
		}
		printf("\n");
	}
}

void splitgraph_GPU(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, \
	hunyuangraph_graph_t **r_lgraph, hunyuangraph_graph_t **r_rgraph, int flag)
{
	int nvtxs, nedges;

	hunyuangraph_graph_t *lgraph,*rgraph;
	lgraph=hunyuangraph_create_cpu_graph();
	rgraph=hunyuangraph_create_cpu_graph();

	int lnvtxs, rnvtxs, lnedges, rnedges;
	int *lmap, *rmap, *ltxadj, *rtxadj;

	nvtxs  = graph->nvtxs;
	nedges = graph->nedges;

	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	// cudaMalloc((void**)&graph->cuda_xadj,sizeof(int) * (graph->nvtxs + 1));
	// cudaMalloc((void**)&graph->cuda_adjncy,sizeof(int) * graph->nedges);
	// cudaMalloc((void**)&graph->cuda_adjwgt,sizeof(int) * graph->nedges);
	// cudaMalloc((void**)&graph->cuda_vwgt,sizeof(int) * graph->nvtxs);
	// cudaMalloc((void**)&graph->cuda_where,sizeof(int) * graph->nvtxs);
	cudaMalloc((void**)&lmap,sizeof(int) * graph->nvtxs);
	cudaMalloc((void**)&rmap,sizeof(int) * graph->nvtxs);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;

	cudaDeviceSynchronize();
	gettimeofday(&begin_memcpy_split, NULL);
	// cudaMemcpy(graph->cuda_xadj,graph->xadj,sizeof(int) * (graph->nvtxs + 1), cudaMemcpyHostToDevice);
	// cudaMemcpy(graph->cuda_adjncy,graph->adjncy,sizeof(int) * graph->nedges, cudaMemcpyHostToDevice);
	// cudaMemcpy(graph->cuda_adjwgt,graph->adjwgt,sizeof(int) * graph->nedges, cudaMemcpyHostToDevice);
	// cudaMemcpy(graph->cuda_vwgt,graph->vwgt,sizeof(int) * graph->nvtxs, cudaMemcpyHostToDevice);
	// cudaMemcpy(graph->cuda_where,graph->where,sizeof(int) * graph->nvtxs, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	gettimeofday(&end_memcpy_split, NULL);
	memcpy_split += (end_memcpy_split.tv_sec - begin_memcpy_split.tv_sec) * 1000 + (end_memcpy_split.tv_usec - begin_memcpy_split.tv_usec) / 1000.0;

	// rgraph:where=0, rgraph:where=1
	// 计算左右子图各自的顶点数
	rnvtxs = thrust::reduce(thrust::device, graph->cuda_where, graph->cuda_where + nvtxs);
	lnvtxs = nvtxs - rnvtxs;

	// lgraph
	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&ltxadj,sizeof(int) * (nvtxs + 1));
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	compute_lnedges<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_where,ltxadj);

	thrust::exclusive_scan(thrust::device, ltxadj, ltxadj + nvtxs + 1, ltxadj);
	cudaMemcpy(&lnedges, &ltxadj[nvtxs],sizeof(int), cudaMemcpyDeviceToHost);

	// printf("lgraph ltxadj\n");

	compute_lmap<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,lmap);
	thrust::inclusive_scan(thrust::device, lmap, lmap + nvtxs, lmap);

	// printf("lgraph map\n");

	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&lgraph->cuda_xadj,sizeof(int) * (lnvtxs + 1));
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	set_lxadj<<<(nvtxs + 128) / 128, 128>>>(nvtxs,lnvtxs,lnedges,graph->cuda_xadj,graph->cuda_where,lmap,ltxadj,lgraph->cuda_xadj);

	// printf("lgraph xadj\n");

	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&lgraph->cuda_adjncy,sizeof(int) * lnedges);
	cudaMalloc((void**)&lgraph->cuda_adjwgt,sizeof(int) * lnedges);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	set_ladjncy_ladjwgt<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
		graph->cuda_where,lgraph->cuda_xadj,lgraph->cuda_adjncy,lgraph->cuda_adjwgt,lmap);
	
	// printf("lgraph adjncy\n");

	// vwgt
	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&lgraph->cuda_vwgt,sizeof(int) * lnvtxs);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	set_lvwgt<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,graph->cuda_vwgt,lmap,lgraph->cuda_vwgt);

	// printf("lgraph vwgt\n");

	// label
	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&lgraph->cuda_label,sizeof(int) * lnvtxs);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	if(flag == 1) set_llabel0<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,lmap,lgraph->cuda_label);
	else 
	{
		cudaMalloc((void**)&graph->cuda_label,sizeof(int) * graph->nvtxs);
		cudaMemcpy(graph->cuda_label,graph->label,sizeof(int) * graph->nvtxs, cudaMemcpyHostToDevice);
		set_llabel1<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,lmap,graph->cuda_label,lgraph->cuda_label);
	}
	// printf("lgraph end\n");

	/*printf("GPU_lgraph:\n");
	cudaDeviceSynchronize();
	exam_graph<<<1,1>>>(lnvtxs,lnedges,lgraph->cuda_xadj,lgraph->cuda_adjncy,lgraph->cuda_vwgt,lgraph->cuda_label);
	cudaDeviceSynchronize();*/

	// rgraph
	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&rtxadj,sizeof(int) * (nvtxs + 1));
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	compute_rnedges<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_where,rtxadj);

	thrust::exclusive_scan(thrust::device, rtxadj, rtxadj + nvtxs + 1, rtxadj);
	cudaMemcpy(&rnedges, &rtxadj[nvtxs],sizeof(int), cudaMemcpyDeviceToHost);

	// printf("rgraph rtxadj\n");

	compute_rmap<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,rmap);
	thrust::inclusive_scan(thrust::device, rmap, rmap + nvtxs, rmap);

	// printf("rgraph map\n");

	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&rgraph->cuda_xadj,sizeof(int) * (rnvtxs + 1));
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	set_rxadj<<<(nvtxs + 128) / 128, 128>>>(nvtxs,rnvtxs,rnedges,graph->cuda_xadj,graph->cuda_where,rmap,rtxadj,rgraph->cuda_xadj);

	// printf("rgraph xadj\n");

	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&rgraph->cuda_adjncy,sizeof(int) * rnedges);
	cudaMalloc((void**)&rgraph->cuda_adjwgt,sizeof(int) * rnedges);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	set_radjncy_radjwgt<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
		graph->cuda_where,rgraph->cuda_xadj,rgraph->cuda_adjncy,rgraph->cuda_adjwgt,rmap);
	
	// printf("rgraph adjncy\n");
	
	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&rgraph->cuda_vwgt,sizeof(int) * rnvtxs);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	set_rvwgt<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,graph->cuda_vwgt,rmap,rgraph->cuda_vwgt);

	// printf("rgraph vwgt\n");

	cudaDeviceSynchronize();
	gettimeofday(&begin_malloc_split, NULL);
	cudaMalloc((void**)&rgraph->cuda_label,sizeof(int) * rnvtxs);
	cudaDeviceSynchronize();
	gettimeofday(&end_malloc_split, NULL);
	malloc_split += (end_malloc_split.tv_sec - begin_malloc_split.tv_sec) * 1000 + (end_malloc_split.tv_usec - begin_malloc_split.tv_usec) / 1000.0;
	if(flag == 1) set_rlabel0<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,rmap,rgraph->cuda_label);
	else set_rlabel1<<<(nvtxs + 127) / 128, 128>>>(nvtxs,graph->cuda_where,rmap,graph->cuda_label,rgraph->cuda_label);

	// printf("rgraph end\n");

  	set_subgraph_tvwgt(graph, lgraph, rgraph);

	/*printf("GPU_rgraph:\n");
	cudaDeviceSynchronize();
	exam_graph<<<1,1>>>(rnvtxs,rnedges,rgraph->cuda_xadj,rgraph->cuda_adjncy,rgraph->cuda_vwgt,rgraph->cuda_label);
	cudaDeviceSynchronize();*/
  	// printf("lgraph->tvwgt=%d lgraph->tvwgt_reverse=%f rgraph->tvwgt=%d rgraph->tvwgt_reverse=%f\n",lgraph->tvwgt[0],lgraph->tvwgt_reverse[0],rgraph->tvwgt[0],rgraph->tvwgt_reverse[0]);
	
	lgraph->xadj = (int *)malloc(sizeof(int) * (lnvtxs + 1));
	lgraph->adjncy = (int *)malloc(sizeof(int) * lnedges);
	lgraph->adjwgt = (int *)malloc(sizeof(int) * lnedges);
	lgraph->vwgt = (int *)malloc(sizeof(int) * lnvtxs);
	lgraph->label = (int *)malloc(sizeof(int) * lnvtxs);

	rgraph->xadj = (int *)malloc(sizeof(int) * (rnvtxs + 1));
	rgraph->adjncy = (int *)malloc(sizeof(int) * rnedges);
	rgraph->adjwgt = (int *)malloc(sizeof(int) * rnedges);
	rgraph->vwgt = (int *)malloc(sizeof(int) * rnvtxs);
	rgraph->label = (int *)malloc(sizeof(int) * rnvtxs);

	cudaDeviceSynchronize();
	gettimeofday(&begin_memcpy_split, NULL);
	cudaMemcpy(lgraph->xadj,lgraph->cuda_xadj,sizeof(int) * (lnvtxs + 1), cudaMemcpyDeviceToHost);
	cudaMemcpy(lgraph->adjncy,lgraph->cuda_adjncy,sizeof(int) * lnedges, cudaMemcpyDeviceToHost);
	cudaMemcpy(lgraph->adjwgt,lgraph->cuda_adjwgt,sizeof(int) * lnedges, cudaMemcpyDeviceToHost);
	cudaMemcpy(lgraph->vwgt,lgraph->cuda_vwgt,sizeof(int) * lnvtxs, cudaMemcpyDeviceToHost);
	cudaMemcpy(lgraph->label,lgraph->cuda_label,sizeof(int) * lnvtxs, cudaMemcpyDeviceToHost);

	cudaMemcpy(rgraph->xadj,rgraph->cuda_xadj,sizeof(int) * (rnvtxs + 1), cudaMemcpyDeviceToHost);
	cudaMemcpy(rgraph->adjncy,rgraph->cuda_adjncy,sizeof(int) * rnedges, cudaMemcpyDeviceToHost);
	cudaMemcpy(rgraph->adjwgt,rgraph->cuda_adjwgt,sizeof(int) * rnedges, cudaMemcpyDeviceToHost);
	cudaMemcpy(rgraph->vwgt,rgraph->cuda_vwgt,sizeof(int) * rnvtxs, cudaMemcpyDeviceToHost);
	cudaMemcpy(rgraph->label,rgraph->cuda_label,sizeof(int) * rnvtxs, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	gettimeofday(&end_memcpy_split, NULL);
	memcpy_split += (end_memcpy_split.tv_sec - begin_memcpy_split.tv_sec) * 1000 + (end_memcpy_split.tv_usec - begin_memcpy_split.tv_usec) / 1000.0;

	cudaDeviceSynchronize();
	gettimeofday(&begin_free_split, NULL);
	cudaFree(lmap);
	cudaFree(rmap);
	cudaFree(ltxadj);
	cudaFree(rtxadj);

	cudaFree(graph->cuda_xadj);
	cudaFree(graph->cuda_adjncy);
	cudaFree(graph->cuda_adjwgt);
	cudaFree(graph->cuda_vwgt);
	cudaFree(graph->cuda_label);
	cudaFree(graph->cuda_where);

	cudaFree(lgraph->cuda_xadj);
	cudaFree(lgraph->cuda_adjncy);
	cudaFree(lgraph->cuda_adjwgt);
	cudaFree(lgraph->cuda_label);
	cudaFree(lgraph->cuda_xadj);

	cudaFree(lgraph->cuda_xadj);
	cudaFree(lgraph->cuda_adjncy);
	cudaFree(lgraph->cuda_adjwgt);
	cudaFree(lgraph->cuda_vwgt);
	cudaFree(lgraph->cuda_label);
	cudaDeviceSynchronize();
	gettimeofday(&end_free_split, NULL);
	free_split += (end_free_split.tv_sec - begin_free_split.tv_sec) * 1000 + (end_free_split.tv_usec - begin_free_split.tv_usec) / 1000.0;

	lgraph->nvtxs  = lnvtxs;
	lgraph->nedges = lnedges;
	rgraph->nvtxs  = rnvtxs;
	rgraph->nedges = rnedges;

	*r_lgraph=lgraph;
	*r_rgraph=rgraph;
}

__global__ void set_part0(int nvtxs, int *part, int *where, int fpart)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		part[ii] = where[ii] + fpart;
	}
}

__global__ void set_part1(int nvtxs, int *part, int *where, int *label, int fpart)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		part[label[ii]] = where[ii] + fpart;
	}
}

/*Cpu Multilevel resursive bisection*/
int hunyuangraph_mlevel_rbbisection(hunyuangraph_admin_t *hunyuangraph_admin, \
	hunyuangraph_graph_t *graph, int nparts, int *part, float *tpwgts, int fpart, int level)
{
	int i,nvtxs,objval;
	int *label,*where;

	hunyuangraph_graph_t *lgraph,*rgraph;
	float wsum,*tpwgts2;

	if(graph->nvtxs==0){
		printf("****You are trying to partition too many parts!****\n");
		return 0;
	}

	nvtxs=graph->nvtxs;

	tpwgts2=hunyuangraph_float_malloc_space(hunyuangraph_admin);
	CPU_malloc(sizeof(float) * 2);
	tpwgts2[0]=hunyuangraph_float_sum((nparts>>1),tpwgts);
	tpwgts2[1]=1.0-tpwgts2[0];

  	objval=hunyuangraph_cpu_mlevelbisect(hunyuangraph_admin,graph,tpwgts2);

	// printf("hunyuangraph_cpu_mlevelbisect\n");

    level++;

	if(level == 1)
	{
		where = graph->where;
		for(i = 0;i < nvtxs;i++)
			part[i] = where[i] + fpart;
	}
	else
	{
		label = graph->label;
		where = graph->where;

		for(i = 0;i < nvtxs;i++){
			part[label[i]] = where[i] + fpart;
		}
	}

	// printf("label\n");

	if(nparts>2){ 
		cudaDeviceSynchronize();
		gettimeofday(&begin_part_slipt,NULL);
		// if(graph->nvtxs > 10000) splitgraph_GPU(hunyuangraph_admin,graph,&lgraph,&rgraph,level);
		// else
		// {
			if(level != 1) hunyuangraph_splitgraph(hunyuangraph_admin,graph,&lgraph,&rgraph);
			else hunyuangraph_splitgraph_first(hunyuangraph_admin,graph,&lgraph,&rgraph);
		// }
		cudaDeviceSynchronize();
		gettimeofday(&end_part_slipt,NULL);
		part_slipt += (end_part_slipt.tv_sec - begin_part_slipt.tv_sec) * 1000 + (end_part_slipt.tv_usec - begin_part_slipt.tv_usec) / 1000.0;
		// printf("hunyuangraph_cpu_mlevelbisect\n");
	}
	
	hunyuangraph_free_graph(&graph);

	wsum=hunyuangraph_float_sum((nparts>>1),tpwgts);
	
	hunyuangraph_tpwgts_rescale((nparts>>1),1.0/wsum,tpwgts);
	hunyuangraph_tpwgts_rescale(nparts-(nparts>>1),1.0/(1.0-wsum),tpwgts+(nparts>>1));

	if(nparts>3){

		objval+=hunyuangraph_mlevel_rbbisection(hunyuangraph_admin,lgraph,(nparts>>1),part,tpwgts,fpart,level);

		objval+=hunyuangraph_mlevel_rbbisection(hunyuangraph_admin,rgraph,nparts-(nparts>>1),part,tpwgts+(nparts>>1),fpart+(nparts>>1),level);
	}
	else if(nparts==3){
		hunyuangraph_free_graph(&lgraph);
		objval+=hunyuangraph_mlevel_rbbisection(hunyuangraph_admin,rgraph,nparts-(nparts>>1),part,tpwgts+(nparts>>1),fpart+(nparts>>1),level);
	}
	
	return objval;

}

/*Set kway balance params*/
void hunyuangraph_set_kway_bal(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
	for(int i = 0, j = 0;i < hunyuangraph_admin->nparts;i++)
		hunyuangraph_admin->part_balance[i + j] = graph->tvwgt_reverse[j] / hunyuangraph_admin->tpwgts[i + j];
}

/*Cpu graph partition algorithm*/
int hunyuangraph_rbbisection(int *nvtxs, int *xadj, int *adjncy, int *vwgt,int *adjwgt, int *nparts, float *tpwgts, float *ubvec, int *objval, int *part, int *tvwgt)
{
  hunyuangraph_graph_t *graph;
  hunyuangraph_admin_t *hunyuangraph_admin;

  hunyuangraph_admin = hunyuangraph_set_graph_admin( *nparts, tpwgts, ubvec);
  CPU_malloc(sizeof(hunyuangraph_admin_t) + sizeof(int) + sizeof(float));

    // cudaDeviceSynchronize();
    // gettimeofday(&begin_set_cpu_graph,NULL);
    graph = hunyuangraph_set_graph(hunyuangraph_admin, *nvtxs, xadj, adjncy, vwgt, adjwgt, tvwgt);
	CPU_malloc(sizeof(hunyuangraph_graph_t) + sizeof(float) * (*nparts * 2 + 1));
    // cudaDeviceSynchronize();
    // gettimeofday(&end_set_cpu_graph,NULL);
    // set_cpu_graph += (end_set_cpu_graph.tv_sec - begin_set_cpu_graph.tv_sec) * 1000 + (end_set_cpu_graph.tv_usec - begin_set_cpu_graph.tv_usec) / 1000.0;

  hunyuangraph_allocatespace(hunyuangraph_admin, graph);           
  
  *objval = hunyuangraph_mlevel_rbbisection(hunyuangraph_admin, graph, *nparts, part, hunyuangraph_admin->tpwgts, 0, 0);
  
  return 1;
 
}

/*Graph initial partition algorithm*/
void hunyuangarph_initialpartition(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
  int objval=0;
  int *bestwhere=NULL;
  float *ubvec=NULL;

  graph->where=(int *)malloc(sizeof(int)*graph->nvtxs);
  CPU_malloc(sizeof(int)*graph->nvtxs);
  hunyuangraph_admin->ncuts = hunyuangraph_admin->nIparts;

  ubvec=(float*)malloc(sizeof(float));
  CPU_malloc(sizeof(float));
  ubvec[0]=(float)pow(hunyuangraph_admin->ubfactors[0],1.0/log(hunyuangraph_admin->nparts));
  
  hunyuangraph_rbbisection(&graph->nvtxs,graph->xadj,graph->adjncy,graph->vwgt,graph->adjwgt, \
    &hunyuangraph_admin->nparts,hunyuangraph_admin->tpwgts,ubvec,&objval,graph->where,graph->tvwgt);
  
  free(ubvec);
  free(bestwhere);
}

/*CUDA-init pwgts array*/
__global__ void initpwgts(int *cuda_pwgts, int nparts)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nparts)
    cuda_pwgts[ii] = 0;
}

/*Compute sum of pwgts*/
__global__ void Sumpwgts(int *cuda_pwgts, int *cuda_where, int *cuda_vwgt, int nvtxs)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nvtxs)
  {
    int t = cuda_where[ii];
    atomicAdd(&cuda_pwgts[t],cuda_vwgt[ii]);
  }
}

__global__ void calculateSum(int nvtxs, int nparts, int *pwgts, int *where, int *vwgt)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	extern __shared__ int cache_d[];
	for (int i = threadIdx.x; i < nparts; i += 128)
		cache_d[i] = 0;
	__syncthreads();

	int t;
	if(ii < nvtxs)
	{
		t = where[ii];
		atomicAdd(&cache_d[t],vwgt[ii]);
	}
	__syncthreads();

	int val;
	for (int i = threadIdx.x; i < nparts; i += 128)
	{
		val = cache_d[i];
		if(val > 0)
		{
			atomicAdd(&pwgts[i],val);
		}
	}
}

/*CUDA-init pwgts array*/
__global__ void inittpwgts(float *tpwgts, float temp, int nparts)
{
	int ii = blockIdx.x*blockDim.x+threadIdx.x;

	if(ii < nparts)
		tpwgts[ii] = temp;
}

__global__ void exam_pwgts(int nparts, float temp, int tvwgt,int *pwgts)
{
	printf("[%f %f]\n",temp / 1.03 * tvwgt, temp * 1.03 * tvwgt);
	for(int i = 0;i < nparts;i++)
		printf("%d ",pwgts[i]);
	printf("\n");
}

/*Malloc initial partition phase to refine phase params*/
void Mallocinit_refineinfo(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
	int nvtxs  = graph->nvtxs;
	int nparts = hunyuangraph_admin->nparts;
	int num    = 0;

	// cudaMalloc((void**)&graph->cuda_where,nvtxs * sizeof(int));
	// cudaMalloc((void**)&graph->cuda_bnd,nvtxs * sizeof(int));
	// cudaMalloc((void**)&graph->cuda_bndnum,sizeof(int));
	// cudaMalloc((void**)&graph->cuda_pwgts,nparts * sizeof(int));
	// cudaMalloc((void**)&graph->cuda_tpwgts,nparts * sizeof(float));
	// cudaMalloc((void**)&graph->cuda_maxwgt,nparts * sizeof(int));
	// cudaMalloc((void**)&graph->cuda_minwgt,nparts * sizeof(int));
	
	graph->cuda_where  = (int *)lmalloc_with_check(sizeof(int) * nvtxs,"where");
	graph->cuda_bnd    = (int *)lmalloc_with_check(sizeof(int) * nvtxs,"bnd");
	graph->cuda_pwgts  = (int *)lmalloc_with_check(sizeof(int) * nparts,"pwgts");
	graph->cuda_tpwgts = (float *)lmalloc_with_check(sizeof(float) * nparts,"tpwgts");
	graph->cuda_maxwgt = (int *)lmalloc_with_check(sizeof(int) * nparts,"maxwgt");
	graph->cuda_minwgt = (int *)lmalloc_with_check(sizeof(int) * nparts,"minwgt");
	graph->cuda_bndnum = (int *)lmalloc_with_check(sizeof(int),"bndnum");

    graph->cuda_bn  = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_bn");
	graph->cuda_bt  = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_bt");
	// graph->cuda_g   = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_g");
	graph->cuda_csr = (int *)lmalloc_with_check(sizeof(int) * 2,"cu_csr");
	graph->cuda_que = (int *)lmalloc_with_check(sizeof(int) * hunyuangraph_admin->nparts * 2,"cu_que");

	cudaMemcpy(graph->cuda_where,graph->where,nvtxs * sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(graph->cuda_bndnum,&num,sizeof(int),cudaMemcpyHostToDevice);

	initpwgts<<<nparts / 32 + 1,32>>>(graph->cuda_pwgts,nparts);

	cudaDeviceSynchronize();
	gettimeofday(&begin_krefine_atomicadd,NULL);
//   Sumpwgts<<<nvtxs/32+1,32>>>(graph->cuda_pwgts,graph->cuda_where,graph->cuda_vwgt,nvtxs);
	calculateSum<<<(nvtxs + 127) / 128,128,nparts * sizeof(int)>>>(nvtxs,nparts,graph->cuda_pwgts,graph->cuda_where,graph->cuda_vwgt);
	cudaDeviceSynchronize();
	gettimeofday(&end_krefine_atomicadd,NULL);
	krefine_atomicadd += (end_krefine_atomicadd.tv_sec - begin_krefine_atomicadd.tv_sec) * 1000 + (end_krefine_atomicadd.tv_usec - begin_krefine_atomicadd.tv_usec) / 1000.0;

	/*if(graph->where == NULL)
	{
		graph->where = (int *)malloc(sizeof(int) * nvtxs);
		cudaMemcpy(graph->where,graph->cuda_where,sizeof(int) * nvtxs,cudaMemcpyDeviceToHost);
	}
	if(graph->vwgt == NULL)
	{
		graph->vwgt = (int *)malloc(sizeof(int) * nvtxs);
		cudaMemcpy(graph->vwgt,graph->cuda_vwgt,sizeof(int) * nvtxs,cudaMemcpyDeviceToHost);
	}
	if(graph->pwgts == NULL)
	{
		graph->pwgts = (int *)malloc(sizeof(int) * nparts);
		for(int i = 0;i < nparts;i++)
			graph->pwgts[i] = 0;
	}

	for(int i = 0;i < nvtxs;i++)
	{
		graph->pwgts[graph->where[i]] += graph->vwgt[i];
	}
	for(int i = 0;i < nparts;i++)
		printf("%d ",graph->pwgts[i]);
	printf("\n");

	cudaDeviceSynchronize();
	exam_where<<<1,1>>>(nparts,graph->cuda_pwgts);
	cudaDeviceSynchronize();*/

	inittpwgts<<<nparts / 32 + 1,32>>>(graph->cuda_tpwgts,hunyuangraph_admin->tpwgts[0],nparts);

	// cudaDeviceSynchronize();
	// exam_pwgts<<<1,1>>>(nparts,hunyuangraph_admin->tpwgts[0],graph->tvwgt[0],graph->cuda_pwgts);
	// cudaDeviceSynchronize();
}

/*Malloc refine params*/
void hunyuangraph_malloc_refineinfo(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
	int nvtxs  = graph->nvtxs;
	int nparts = hunyuangraph_admin->nparts;
	int num    = 0;

	// cudaMalloc((void**)&graph->cuda_bnd,nvtxs * sizeof(int));
	// cudaMalloc((void**)&graph->cuda_bndnum,sizeof(int));
	// cudaMalloc((void**)&graph->cuda_pwgts,nparts * sizeof(int));
	// cudaMalloc((void**)&graph->cuda_tpwgts,nparts * sizeof(float));
	// cudaMalloc((void**)&graph->cuda_maxwgt,nparts * sizeof(int));
	// cudaMalloc((void**)&graph->cuda_minwgt,nparts * sizeof(int));

	graph->cuda_bnd    = (int *)lmalloc_with_check(sizeof(int) * nvtxs,"bnd");
	graph->cuda_pwgts  = (int *)lmalloc_with_check(sizeof(int) * nparts,"pwgts");
	graph->cuda_tpwgts = (float *)lmalloc_with_check(sizeof(float) * nparts,"tpwgts");
	graph->cuda_maxwgt = (int *)lmalloc_with_check(sizeof(int) * nparts,"maxwgt");
	graph->cuda_minwgt = (int *)lmalloc_with_check(sizeof(int) * nparts,"minwgt");
	graph->cuda_bndnum = (int *)lmalloc_with_check(sizeof(int),"bndnum");

	graph->cuda_bn  = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_bn");
	graph->cuda_bt  = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_bt");
	// graph->cuda_g   = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_g");
	graph->cuda_csr = (int *)lmalloc_with_check(sizeof(int) * 2,"cu_csr");
	graph->cuda_que = (int *)lmalloc_with_check(sizeof(int) * hunyuangraph_admin->nparts * 2,"cu_que");

	cudaMemcpy(graph->cuda_bndnum,&num,sizeof(int),cudaMemcpyHostToDevice);

	initpwgts<<<nparts/32+1,32>>>(graph->cuda_pwgts,nparts);

	cudaDeviceSynchronize();
	gettimeofday(&begin_krefine_atomicadd,NULL);
	// Sumpwgts<<<nvtxs/32+1,32>>>(graph->cuda_pwgts,graph->cuda_where,graph->cuda_vwgt,nvtxs);
	calculateSum<<<(nvtxs + 127) / 128,128,nparts * sizeof(int)>>>(nvtxs,nparts,graph->cuda_pwgts,graph->cuda_where,graph->cuda_vwgt);
	cudaDeviceSynchronize();
	gettimeofday(&end_krefine_atomicadd,NULL);
	krefine_atomicadd += (end_krefine_atomicadd.tv_sec - begin_krefine_atomicadd.tv_sec) * 1000 + (end_krefine_atomicadd.tv_usec - begin_krefine_atomicadd.tv_usec) / 1000.0;

	/*if(graph->where == NULL)
	{
		graph->where = (int *)malloc(sizeof(int) * nvtxs);
		cudaMemcpy(graph->where,graph->cuda_where,sizeof(int) * nvtxs,cudaMemcpyDeviceToHost);
	}
	if(graph->vwgt == NULL)
	{
		graph->vwgt = (int *)malloc(sizeof(int) * nvtxs);
		cudaMemcpy(graph->vwgt,graph->cuda_vwgt,sizeof(int) * nvtxs,cudaMemcpyDeviceToHost);
	}
	if(graph->pwgts == NULL)
	{
		graph->pwgts = (int *)malloc(sizeof(int) * nparts);
		for(int i = 0;i < nparts;i++)
			graph->pwgts[i] = 0;
	}

	for(int i = 0;i < nvtxs;i++)
	{
		graph->pwgts[graph->where[i]] += graph->vwgt[i];
	}
	for(int i = 0;i < nparts;i++)
		printf("%d ",graph->pwgts[i]);
	printf("\n");

	cudaDeviceSynchronize();
	exam_where<<<1,1>>>(nparts,graph->cuda_pwgts);
	cudaDeviceSynchronize();*/

	inittpwgts<<<nparts/32+1,32>>>(graph->cuda_tpwgts,hunyuangraph_admin->tpwgts[0],nparts);

	// cudaDeviceSynchronize();
	// exam_pwgts<<<1,1>>>(nparts,hunyuangraph_admin->tpwgts[0],graph->tvwgt[0],graph->cuda_pwgts);
	// cudaDeviceSynchronize();
}

/*CUDA-get the max/min pwgts*/
__global__ void Sum_maxmin_pwgts(int *maxwgt, int *minwgt, float *tpwgts, int tvwgt, int nparts)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nparts)
	{
		float result = tpwgts[ii] * tvwgt;

		maxwgt[ii] = int(result * 1.03);
		minwgt[ii] = int(result / 1.03);
	}
}

/*CUDA-init boundary vertex num*/
__global__ void initbndnum(int *bndnum)
{
  bndnum[0] = 0;
}

/*CUDA-find vertex where ed-id>0 */
__global__ void Find_real_bnd_info(int *cuda_real_bnd_num, int *cuda_real_bnd, int *where, \
	int *xadj, int *adjncy, int *adjwgt, int nvtxs, int nparts)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs) // && moved[ii] == 0
	{
		int i, k, begin, end, me, other, from;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		me    = 0;
		other = 0;
		from  = where[ii];

		for(i = begin;i < end;i++)
		{
			k = adjncy[i];
			if(where[k] == from) me += adjwgt[i];
			else other += adjwgt[i];
		}
		if(other > me) cuda_real_bnd[atomicAdd(&cuda_real_bnd_num[0],1)] = ii;
	}
}

/**/
__global__ void init_bnd_info(int *bnd_info, int length)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < length)
    bnd_info[ii] = 0;
}

/*CUDA-find boundary vertex should ro which part*/
__global__ void find_kayparams(int *cuda_real_bnd_num, int *bnd_info, int *cuda_real_bnd, int *where, \
int *xadj, int *adjncy, int *adjwgt, int nparts, int *cuda_bn, int *cuda_bt)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii<cuda_real_bnd_num[0])
  {
    int pi, other, i, k, me_wgt, other_wgt;
    int start, end, begin, last;

    pi    = cuda_real_bnd[ii];
    begin = xadj[pi];
    last  = xadj[pi+1];
    start = nparts * ii;
    end   = nparts * (ii+1);
    other = where[pi];

    for(i = begin;i < last;i++)
    {
      k = adjncy[i];
      k = start + where[k];
      bnd_info[k] += adjwgt[i];
    }

    me_wgt = other_wgt = bnd_info[start + other];

    for(i=start;i<end;i++)
    {
      k = bnd_info[i];
      if(k > other_wgt)
      {
        other_wgt = k;
        other     = i - start;
      }
    }

    // cuda_g[ii]  = other_wgt - me_wgt;
    cuda_bt[ii] = other;
    cuda_bn[ii] = pi;

	// atomicAdd(&cuda_que[other],1);

  }
}

/*CUDA-init params*/
__global__ void initcucsr(int *cu_csr, int *bndnum)
{
  cu_csr[0] = 0;
  cu_csr[1] = bndnum[0];
}

/**/
__global__ void init_cu_que(int *cuda_que, int length)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < length)
    cuda_que[ii] = 0;
}

/*CUDA-get a csr array*/
__global__ void findcsr(int *cuda_bt, int *cuda_que, int *bnd_num, int nparts)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nparts)
  {
    int i, t;
    int begin, end;
    begin = 2 * ii;
    end   = bnd_num[0];

    for(i = 0;i < end;i++)
    {
      if(ii == cuda_bt[i])
      {
        cuda_que[begin] = i;
        break; 
      }
    }

    t = cuda_que[begin];

    if(t!=-1)
    {
      for(i = t;i < end;i++)
      {
        if(cuda_bt[i] != ii)
        {
          cuda_que[begin + 1] = i - 1;
          break; 
        }
      }
    }

    t = 2 * cuda_bt[end - 1] + 1;
    cuda_que[t] = end - 1;
  }
}

__global__ void findcsr_me(int *cuda_bt, int *cuda_que, int *bnd_num)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < bnd_num[0])
	{
		int to = cuda_bt[ii];
		atomicAdd(&cuda_que[to],1);
	}
}

__global__ void calculateCsr(int *bnd_num, int nparts, int *cuda_que, int *cuda_bt)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	extern __shared__ int cache_d[];
	for (int i = threadIdx.x; i < nparts; i += 128)
		cache_d[i] = 0;
	__syncthreads();

	int t;
	if(ii < bnd_num[0])
	{
		t = cuda_bt[ii];
		atomicAdd(&cache_d[t],1);
	}
	__syncthreads();

	int val;
	for (int i = threadIdx.x; i < nparts; i += 128)
	{
		val = cache_d[i];
		if(val > 0)
		{
			atomicAdd(&cuda_que[i],val);
		}
	}
}

__global__ void exam_que(int nparts, int *cuda_que)
{
	for(int i = 0;i <= nparts;i++)
		printf("i=%5d que=%5d\n",i,cuda_que[i]);
}

/*Find boundary vertex information*/
void hunyuangraph_findgraphbndinfo(hunyuangraph_admin_t *hunyuangraph_admin,hunyuangraph_graph_t *graph)
{
	int nvtxs  = graph->nvtxs;
	int nparts = hunyuangraph_admin->nparts;
	int bnd_num; 

	initbndnum<<<1,1>>>(graph->cuda_bndnum);

	cudaDeviceSynchronize();
	gettimeofday(&begin_bndinfo,NULL);
	Find_real_bnd_info<<<nvtxs / 32 + 1,32>>>(graph->cuda_bndnum,graph->cuda_bnd,graph->cuda_where,\
		graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,nvtxs,nparts);
	cudaDeviceSynchronize();
	gettimeofday(&end_bndinfo,NULL);
	bndinfo_Find_real_bnd_info += (end_bndinfo.tv_sec - begin_bndinfo.tv_sec) * 1000 + (end_bndinfo.tv_usec - begin_bndinfo.tv_usec) / 1000.0; 
  
	cudaMemcpy(&bnd_num,graph->cuda_bndnum, sizeof(int), cudaMemcpyDeviceToHost);
  
	if(bnd_num > 0)
	{
		// cudaMalloc((void**)&graph->cuda_info, bnd_num * nparts * sizeof(int));
		graph->cuda_info = (int *)rmalloc_with_check(sizeof(int) * bnd_num * nparts,"info");

		cudaDeviceSynchronize();
		gettimeofday(&begin_bndinfo,NULL);
		init_bnd_info<<<bnd_num * nparts / 32 + 1,32>>>(graph->cuda_info, bnd_num * nparts);
		cudaDeviceSynchronize();
		gettimeofday(&end_bndinfo,NULL);
		bndinfo_init_bnd_info += (end_bndinfo.tv_sec - begin_bndinfo.tv_sec) * 1000 + (end_bndinfo.tv_usec - begin_bndinfo.tv_usec) / 1000.0; 

		cudaDeviceSynchronize();
		gettimeofday(&begin_bndinfo,NULL);
		find_kayparams<<<bnd_num/32+1,32>>>(graph->cuda_bndnum,graph->cuda_info,graph->cuda_bnd,graph->cuda_where,\
			graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,nparts,graph->cuda_bn,graph->cuda_bt);
		cudaDeviceSynchronize();
		gettimeofday(&end_bndinfo,NULL);
		bndinfo_find_kayparams += (end_bndinfo.tv_sec - begin_bndinfo.tv_sec) * 1000 + (end_bndinfo.tv_usec - begin_bndinfo.tv_usec) / 1000.0; 

		cudaDeviceSynchronize();
		gettimeofday(&begin_bndinfo,NULL);
		initcucsr<<<1,1>>>(graph->cuda_csr,graph->cuda_bndnum);
		cudaDeviceSynchronize();
		gettimeofday(&end_bndinfo,NULL);
		bndinfo_initcucsr += (end_bndinfo.tv_sec - begin_bndinfo.tv_sec) * 1000 + (end_bndinfo.tv_usec - begin_bndinfo.tv_usec) / 1000.0; 

		cudaDeviceSynchronize();
		gettimeofday(&begin_bndinfo,NULL);
		bb_segsort(graph->cuda_bt, graph->cuda_bn, bnd_num, graph->cuda_csr, 1);
		// thrust::sort(thrust::device, kp, kp + nvtxs, compRule());
		cudaDeviceSynchronize();
		gettimeofday(&end_bndinfo,NULL);
		bndinfo_bb_segsort += (end_bndinfo.tv_sec - begin_bndinfo.tv_sec) * 1000 + (end_bndinfo.tv_usec - begin_bndinfo.tv_usec) / 1000.0; 

		cudaDeviceSynchronize();
		gettimeofday(&begin_bndinfo,NULL);
		init_cu_que<<<2 * nparts / 32 + 1,32>>>(graph->cuda_que, 2 * nparts);
		// init_cu_que<<<(nparts + 32) / 32,32>>>(graph->cuda_que, nparts + 1);
		cudaDeviceSynchronize();
		gettimeofday(&end_bndinfo,NULL);
		bndinfo_init_cu_que += (end_bndinfo.tv_sec - begin_bndinfo.tv_sec) * 1000 + (end_bndinfo.tv_usec - begin_bndinfo.tv_usec) / 1000.0; 

		cudaDeviceSynchronize();
		gettimeofday(&begin_bndinfo,NULL);
		// findcsr<<<nparts/32+1,32>>>(graph->cuda_bt,graph->cuda_que,graph->cuda_bndnum,nparts);

		// findcsr_me<<<bnd_num/32+1,32>>>(graph->cuda_bt,graph->cuda_que,graph->cuda_bndnum);

		// init_cu_que<<<(nparts + 31) / 32,32>>>(graph->cuda_que, nparts + 1);
		calculateCsr<<<(bnd_num + 127) / 128,128,nparts * sizeof(int)>>>(graph->cuda_bndnum,nparts,graph->cuda_que,graph->cuda_bt);
		
		// printf("shared:\n");
		// cudaDeviceSynchronize();
		// exam_que<<<1,1>>>(nparts,graph->cuda_que);
		// cudaDeviceSynchronize();

		thrust::exclusive_scan(thrust::device, graph->cuda_que, graph->cuda_que + nparts + 1, graph->cuda_que);
		cudaDeviceSynchronize();
		gettimeofday(&end_bndinfo,NULL);
		bndinfo_findcsr += (end_bndinfo.tv_sec - begin_bndinfo.tv_sec) * 1000 + (end_bndinfo.tv_usec - begin_bndinfo.tv_usec) / 1000.0; 

		rfree_with_check(sizeof(int) * bnd_num * nparts,"info");
		// cudaFree(graph->cuda_info);
	}

	graph->cpu_bndnum=(int *)malloc(sizeof(int));
	graph->cpu_bndnum[0]=bnd_num;
}

/*CUDA-move vertex*/
__global__ void Exnode_part1(int *cuda_que, int *pwgts, int *bnd, int *bndto, int *vwgt,\
  int *maxvwgt, int *minvwgt, int *where)
{
	int ii  = blockIdx.x;
	int iii = threadIdx.x;

	if(iii == 0)
	{
		int i, me, to, vvwgt;
		int pvmax, pvmin, mepwgts, topwgts;
		int begin, end, k;
		
		begin = cuda_que[ii];
		end   = cuda_que[ii + 1];

		pvmax   = maxvwgt[me];
		pvmin   = minvwgt[me];

		for(i = begin;i < end;i++)
		{
			k     = bnd[i];
			me    = where[k];
			to    = bndto[i];

			if(me < to)
			{
				vvwgt   = vwgt[k];
				mepwgts = pwgts[me];
				topwgts = pwgts[to];

				if(((topwgts + vvwgt >= pvmin) && (topwgts + vvwgt <= pvmax))\
					&&((mepwgts - vvwgt >= pvmin) && (mepwgts - vvwgt <= pvmax)))
				{
					atomicAdd(&pwgts[to],vvwgt);
					atomicSub(&pwgts[me],vvwgt);
					where[k] = to;
				}
			}
		}
	}
}

/*CUDA-move vertex*/
__global__ void Exnode_part2(int *cuda_que, int *pwgts, int *bnd, int *bndto, int *vwgt,\
  int *maxvwgt, int *minvwgt, int *where)
{
	int ii  = blockIdx.x;
	int iii = threadIdx.x;

	if(iii == 0)
	{
		int i, me, to, vvwgt;
		int pvmax, pvmin, mepwgts, topwgts;
		int begin, end, k;
		
		begin = cuda_que[ii];
		end   = cuda_que[ii + 1];

		pvmax   = maxvwgt[me];
		pvmin   = minvwgt[me];

		for(i = begin;i < end;i++)
		{
			k     = bnd[i];
			me    = where[k];
			to    = bndto[i];

			if(me > to)
			{
				vvwgt   = vwgt[k];
				mepwgts = pwgts[me];
				topwgts = pwgts[to];

				if(((topwgts + vvwgt >= pvmin) && (topwgts + vvwgt <= pvmax))\
					&&((mepwgts - vvwgt >= pvmin) && (mepwgts - vvwgt <= pvmax)))
				{
					atomicAdd(&pwgts[to],vvwgt);
					atomicSub(&pwgts[me],vvwgt);
					where[k] = to;
				}
			}
		}
	}
}

/*Graph multilevel uncoarsening algorithm*/
void hunyuangraph_k_refinement(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
	int nparts = hunyuangraph_admin->nparts;
	int nvtxs  = graph->nvtxs;

	cudaDeviceSynchronize();
	gettimeofday(&begin_general,NULL);
	Sum_maxmin_pwgts<<<nparts / 32 + 1,32>>>(graph->cuda_maxwgt,graph->cuda_minwgt,graph->cuda_tpwgts,graph->tvwgt[0],nparts);
	cudaDeviceSynchronize();
	gettimeofday(&end_general,NULL);
	uncoarsen_Sum_maxmin_pwgts += (end_general.tv_sec - begin_general.tv_sec) * 1000 + (end_general.tv_usec - begin_general.tv_usec) / 1000.0;

	int i;
	for(i = 0;i < 5;i++)
	{
		hunyuangraph_findgraphbndinfo(hunyuangraph_admin,graph);
		if(graph->cpu_bndnum[0] > 0)
		{
			cudaDeviceSynchronize();
			gettimeofday(&begin_general,NULL);
			Exnode_part1<<<nparts,1>>>(graph->cuda_que,graph->cuda_pwgts,graph->cuda_bn,graph->cuda_bt,graph->cuda_vwgt,\
				graph->cuda_maxwgt,graph->cuda_minwgt,graph->cuda_where);
			cudaDeviceSynchronize();
			gettimeofday(&end_general,NULL);
			uncoarsen_Exnode_part1 += (end_general.tv_sec - begin_general.tv_sec) * 1000 + (end_general.tv_usec - begin_general.tv_usec) / 1000.0;

			cudaDeviceSynchronize();
			gettimeofday(&begin_general,NULL);
			Exnode_part2<<<nparts,1>>>(graph->cuda_que,graph->cuda_pwgts,graph->cuda_bn,graph->cuda_bt,graph->cuda_vwgt,\
				graph->cuda_maxwgt,graph->cuda_minwgt,graph->cuda_where);
			cudaDeviceSynchronize();
			gettimeofday(&end_general,NULL);
			uncoarsen_Exnode_part2 += (end_general.tv_sec - begin_general.tv_sec) * 1000 + (end_general.tv_usec - begin_general.tv_usec) / 1000.0;
		}
		else
			break;
	}
}

__global__ void init_0(int nvtxs, int *num)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
		num[ii] = 0;
}

__global__ void set_bnd(int nvtxs, int *xadj, int *adjncy, int *where, int *moved, int *bnd)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		if(moved[ii] == 0)
		{
			int begin, end, me, flag, i, j;
			begin = xadj[ii];
			end   = xadj[ii + 1];
			me    = where[ii];
			flag  = 0;

			for(i = begin;i < end;i++)
			{
				j  = adjncy[i];
				if(where[j] != me)
				{
					flag = 1;
					break;
				}
			}

			if(flag == 1) bnd[ii] = 1;
			else bnd[ii] = 0;
		}
		else bnd[ii] = 0;
	}
}

__global__ void calculate_to(int nvtxs, int nparts, int *xadj, int *adjncy, int *adjwgt, int *where, int *bnd, int *moveto, int *gain)
{
	int ii  = blockIdx.x;
	int tid = threadIdx.x;

	if(bnd[ii] == 1)
	{
		extern __shared__ int cache_all[];
		int *cache_d   = cache_all;
		int *cache_ptr = cache_all + nparts;
		for (int i = tid; i < nparts; i += 32)
		{
			cache_d[i]   = 0;
			cache_ptr[i] = i;
		}
		__syncthreads();

		int begin, end, me, flag, i, j, id;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		me    = where[ii];
		flag  = 0;

		for(i = begin + tid;i < end;i += 32)
		{
			j  = adjncy[i];
			atomicAdd(&cache_d[where[j]],adjwgt[i]);
		}
		__syncthreads();

		id = cache_d[me];
		__syncthreads();

		// Warp内reduce选择最大度
			// 数据缩小至1个warp
		if(nparts > 32)
		{
			for(i = tid;i < nparts;i += 32)
			{
				if(i + 32 < nparts && cache_d[tid] < cache_d[i + 32])
				{
					cache_d[tid]   = cache_d[i + 32];
					cache_ptr[tid] = cache_ptr[tid + 32];
				}
			}
			__syncthreads();
		}

		if(tid + 16 < nparts && cache_d[tid] < cache_d[tid + 16])
		{
			cache_d[tid]   = cache_d[tid + 16];
			cache_ptr[tid] = cache_ptr[tid + 16];
		}
		__syncthreads();
		if(tid + 8 < nparts && cache_d[tid] < cache_d[tid + 8])
		{
			cache_d[tid]   = cache_d[tid + 8];
			cache_ptr[tid] = cache_ptr[tid + 8];
		}
		__syncthreads();
		if(tid + 4 < nparts && cache_d[tid] < cache_d[tid + 4])
		{
			cache_d[tid]   = cache_d[tid + 4];
			cache_ptr[tid] = cache_ptr[tid + 4];
		}
		__syncthreads();
		if(tid + 2 < nparts && cache_d[tid] < cache_d[tid + 2])
		{
			cache_d[tid]   = cache_d[tid + 2];
			cache_ptr[tid] = cache_ptr[tid + 2];
		}
		__syncthreads();
		if(tid + 1 < nparts && cache_d[tid] < cache_d[tid + 1])
		{
			cache_d[tid]   = cache_d[tid + 1];
			cache_ptr[tid] = cache_ptr[tid + 1];
		}

		if(tid == 0)
		{
			// printf("ii=%d cache=%d ptr=%d where=%d\n",ii,cache_d[0],cache_ptr[0],me);
			flag = cache_ptr[0];
			if(flag == me) moveto[ii] = -1;
			else
			{
				moveto[ii] = flag;
				gain[ii]   = cache_d[0] - id;
			}
		}
	}
	else
	{
		if(tid == 0)
			moveto[ii] = -1;
	}
}

__global__ void init_max(int nparts,int *max, int *maxptr)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nparts)
	{
		max[ii] = 0;
		maxptr[ii] = -1;
	}
}

__global__ void calculate_Max(int nvtxs, int nparts, int *max, int *maxptr, int *moveto, int *gain, int *moved)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;
	
	/*extern __shared__ int cache_d[];
	extern __shared__ int cache_ptr[];
	for (int i = threadIdx.x; i < nparts; i += 128)
	{
		cache_d[i]   = 0;
		cache_ptr[i] = -1;
	}
	__syncthreads();

	if(ii < nvtxs && moved[ii] == -1)
	{
		int t = moveto[ii];
		int val = gain[ii];
		if(t != -1)
		{
			atomicMax(&cache_d[t],val);
			__syncthreads();
			if(cache_d[t] == val) cache_ptr[t] = ii;
		}
		__syncthreads();

		printf("1\n");

			for (int i = 0; i < nparts; i++)
			{
				printf("ii=%d cache=%d cache_ptr=%d\n",ii,cache_d[i],cache_ptr[i]);
			}

		int ptr;
		for (int i = threadIdx.x; i < nparts; i += 128)
		{
			ptr = cache_ptr[i];
			if(ptr != -1)
			{
				atomicMax(&max[t],val);
				if(max[t] == val) maxptr[t] = ptr;
			}
		}
	}*/

	if(ii < nvtxs)
	{
		int t = moveto[ii];
		if(t != -1)
		{
			int val = gain[ii];
			if(val > 0)
			{
				atomicMax(&max[t],val);
				__syncthreads();
				if(max[t] == val) atomicExch(&maxptr[t],ii);
			}
		}
	}
}

__global__ void calculate_move(int nparts, int *maxptr, int *moveto, int *where, int *moved, int *pwgts, int *vwgt)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(ii < nparts)
	{
		int ptr    = maxptr[ii];
		if(ptr != -1)
		{
			int val = vwgt[ptr];
			atomicSub(&pwgts[where[ptr]],val);
			atomicAdd(&pwgts[ii],val);
			where[ptr] = ii;
			moved[ptr] = 1;
		}
	}
}

__global__ void calculate_overweight(int nparts, float temp, int tvwgt,int *pwgts, int *overwgt, int *overthin)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nparts)
	{
		int val;
		float max, min;
		min = temp / 1.03 * tvwgt;
		max = temp * 1.03 * tvwgt;
		val = pwgts[ii];

		if(val > max) overwgt[ii] = 1;
		else overwgt[ii] = 0;

		if(val < min) overthin[ii] = 1;
		else overthin[ii] = 0;
	}
}

__global__ void set_bnd_reba(int nvtxs, int *xadj, int *adjncy, int *where, int *overwgt, int* overthin, int *bnd)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;

	if(ii < nvtxs)
	{
		int me = where[ii];
		if(overwgt[me] == 1)
		{
			int begin, end, flag, i, j, wj;
			begin = xadj[ii];
			end   = xadj[ii + 1];

			for(i = begin;i < end;i++)
			{
				j  = adjncy[i];
				wj = where[j];
				if(wj != me && overthin[wj] == 1)
				{
					flag = 1;
					break;
				}
			}

			if(flag == 1) bnd[ii] = 1;
			else bnd[ii] = 0;
		}
		else bnd[ii] = 0;
	}
}

__global__ void calculate_to_reba(int nvtxs, int nparts, int *xadj, int *adjncy, int *adjwgt, int *where, int *bnd, int *overthin, int *overwgt, int *moveto, int *gain)
{
	int ii  = blockIdx.x;
	int tid = threadIdx.x;

	if(bnd[ii] == 1)
	{
		extern __shared__ int cache_all[];
		int *cache_d   = cache_all;
		int *cache_ptr = cache_all + nparts;
		for (int i = tid; i < nparts; i += 32)
		{
			cache_d[i]   = 0;
			cache_ptr[i] = i;
		}
		__syncthreads();

		int begin, end, me, flag, i, j, wj;
		begin = xadj[ii];
		end   = xadj[ii + 1];
		me    = where[ii];
		flag  = 0;

		for(i = begin + tid;i < end;i += 32)
		{
			j  = adjncy[i];
			wj = where[j];
			if(wj != me && overthin[wj] == 1) atomicAdd(&cache_d[wj],adjwgt[i]);
		}
		__syncthreads();


		/*if(tid == 0)
		{
			for(i = 0;i < nparts;i++)
			{
				printf("ii=%d nparts=%d cache=%d \n",ii,i,cache_d[i]);
			}
		}
		__syncthreads();*/

		// Warp内reduce选择最大度
			// 数据缩小至1个warp
		if(nparts > 32)
		{
			for(i = tid;i < nparts;i += 32)
			{
				if(i + 32 < nparts && cache_d[tid] < cache_d[i + 32])
				{
					cache_d[tid]   = cache_d[i + 32];
					cache_ptr[tid] = cache_ptr[tid + 32];
				}
			}
			__syncthreads();
		}

		if(tid + 16 < nparts && cache_d[tid] < cache_d[tid + 16])
		{
			cache_d[tid]   = cache_d[tid + 16];
			cache_ptr[tid] = cache_ptr[tid + 16];
		}
		__syncthreads();
		if(tid + 8 < nparts && cache_d[tid] < cache_d[tid + 8])
		{
			cache_d[tid]   = cache_d[tid + 8];
			cache_ptr[tid] = cache_ptr[tid + 8];
		}
		__syncthreads();
		if(tid + 4 < nparts && cache_d[tid] < cache_d[tid + 4])
		{
			cache_d[tid]   = cache_d[tid + 4];
			cache_ptr[tid] = cache_ptr[tid + 4];
		}
		__syncthreads();
		if(tid + 2 < nparts && cache_d[tid] < cache_d[tid + 2])
		{
			cache_d[tid]   = cache_d[tid + 2];
			cache_ptr[tid] = cache_ptr[tid + 2];
		}
		__syncthreads();
		if(tid + 1 < nparts && cache_d[tid] < cache_d[tid + 1])
		{
			cache_d[tid]   = cache_d[tid + 1];
			cache_ptr[tid] = cache_ptr[tid + 1];
		}

		if(tid == 0)
		{
			// printf("ii=%d cache=%d ptr=%d where=%d\n",ii,cache_d[0],cache_ptr[0],me);
			flag = cache_ptr[0];
			if(overthin[flag] == 0) moveto[ii] = -1;
			else
			{
				moveto[ii] = cache_ptr[0];
				gain[ii]   = cache_d[0];
			}
		}
	}
	else
	{
		if(tid == 0)
			moveto[ii] = -1;
	}
}

__global__ void calculate_Max_reba(int nvtxs, int nparts, int *max, int *maxptr, int *moveto, int *gain)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(ii < nvtxs)
	{
		int t = moveto[ii];
		if(t != -1)
		{
			int val = gain[ii];
			atomicMax(&max[t],val);
			__syncthreads();
			if(max[t] == val) atomicExch(&maxptr[t],ii);
		}
	}
}

__global__ void calculate_move_reba(int nparts, int *maxptr, int *moveto, int *where, int *pwgts, int *vwgt)
{
	int ii = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(ii < nparts)
	{
		int ptr    = maxptr[ii];
		if(ptr != -1)
		{
			int val = vwgt[ptr];
			atomicSub(&pwgts[where[ptr]],val);
			atomicAdd(&pwgts[ii],val);
			where[ptr] = ii;
		}
	}
}

__global__ void exam_csrpart(int nvtxs, int *xadj, int *adjncy, int *adjwgt, int *where)
{
	for(int i = 0;i < nvtxs;i++)
	{
		for(int j = xadj[i];j < xadj[i + 1];j++)
		{
			printf("%d ",adjncy[j]);
		}
		printf("\n");
		for(int j = xadj[i];j < xadj[i + 1];j++)
		{
			printf("%d ",adjwgt[j]);
		}
		printf("\n");
		for(int j = xadj[i];j < xadj[i + 1];j++)
		{
			printf("%d ",where[adjncy[j]]);
		}
		printf("\n%d ii=%d\n",where[i],i);
		printf("\n");
	}
}

__global__ void exam_gain(int nvtxs,int *where,int *moveto, int *gain)
{
	for(int i = 0;i < nvtxs;i++)
	{
		if(moveto[i] != -1) printf("i=%d from=%d gain[i]=%d moveto=%d\n",i,where[i],gain[i],moveto[i]);
	}
}

__global__ void exam_max(int nparts,int *cu_max, int *cu_maxptr)
{
	for(int i = 0;i < nparts;i++)
	{
		if(cu_maxptr[i] != -1) printf("nparts=%d cu_max=%d cu_maxptr=%d\n",i,cu_max[i],cu_maxptr[i]);
	}
}

__global__ void exam_over(int nparts, int *cu_overwgt, int *cu_overthin)
{
	for(int i = 0;i < nparts;i++)
	{
		printf("nparts=%d cu_overwgt=%d cu_overthin=%d\n",i,cu_overwgt[i],cu_overthin[i]);
	}
}

void hunyuangraph_k_refinement_me(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{
	int nvtxs, nedges, nparts, threshold;
	nvtxs  = graph->nvtxs;
	nedges = graph->nedges;
	nparts = hunyuangraph_admin->nparts;
	threshold = hunyuangraph_max(nvtxs / (6000 * (hunyuangraph_compute_log2(nparts))),8 * nparts);
	printf("threshold=%d\n",threshold);
	
	int *cu_bnd, *cu_moved, *cu_moveto, *cu_gain, *cu_max, *cu_maxptr, *cu_overwgt, *cu_overthin;

	cudaMalloc((void**)&cu_bnd,sizeof(int) * nvtxs);
	cudaMalloc((void**)&cu_moved,sizeof(int) * nvtxs);
	cudaMalloc((void**)&cu_moveto,sizeof(int) * nvtxs);
	cudaMalloc((void**)&cu_gain,sizeof(int) * nvtxs);
	cudaMalloc((void**)&cu_max,sizeof(int) * nparts);
	cudaMalloc((void**)&cu_maxptr,sizeof(int) * nparts);
	cudaMalloc((void**)&cu_overwgt,sizeof(int) * nparts);
	cudaMalloc((void**)&cu_overthin,sizeof(int) * nparts);

	// cudaDeviceSynchronize();
	// exam_csrpart<<<1,1>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,graph->cuda_where);
	// cudaDeviceSynchronize();

	// init_0<<<(nvtxs + 127) / 128,128>>>(nvtxs,cu_bnd);

	// 初始化移动记录数组
	init_0<<<(nvtxs + 127) / 128,128>>>(nvtxs,cu_moved);

	// 先移动
	int iter;
	for(iter = 0;iter < threshold;iter++)
	{
		// printf("iter=%d\n",iter);
		// 标记边界顶点
		set_bnd<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_where,cu_moved,cu_bnd);

		// 计算顶点移动
		calculate_to<<<nvtxs,32,2 * nparts * sizeof(int)>>>(nvtxs,nparts,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,graph->cuda_where,cu_bnd,cu_moveto,cu_gain);

		// cudaDeviceSynchronize();
		// exam_gain<<<1,1>>>(nvtxs,graph->cuda_where,cu_moveto,cu_gain);
		// cudaDeviceSynchronize();

		init_max<<<(nparts + 31) / 32, 32>>>(nparts,cu_max,cu_maxptr);

		// 选择移动至相同分区的最佳顶点
		calculate_Max<<<(nvtxs + 127) / 128,128>>>(nvtxs,nparts,cu_max,cu_maxptr,cu_moveto,cu_gain,cu_moved);

		// cudaDeviceSynchronize();
		// exam_max<<<1,1>>>(nparts,cu_max,cu_maxptr);
		// cudaDeviceSynchronize();

		// 移动并记录至移动记录数组
		calculate_move<<<(nparts + 31) / 32,32>>>(nparts,cu_maxptr,cu_moveto,graph->cuda_where,cu_moved,graph->cuda_pwgts,graph->cuda_vwgt);
	}

	// cudaDeviceSynchronize();
	// exam_pwgts<<<1,1>>>(nparts,hunyuangraph_admin->tpwgts[0],graph->tvwgt[0],graph->cuda_pwgts);
	// cudaDeviceSynchronize();

	// 再平衡
	printf("rbalance\n");
	// 记录过轻或过重 (0:正常,1:过重/过轻)
	calculate_overweight<<<(nparts + 31) / 32, 32>>>(nparts,hunyuangraph_admin->tpwgts[0],graph->tvwgt[0],graph->cuda_pwgts,cu_overwgt,cu_overthin);
	int val_weight = thrust::reduce(thrust::device, cu_overwgt, cu_overwgt + nparts);
	int val_thin   = thrust::reduce(thrust::device, cu_overthin, cu_overthin + nparts);

	// cudaDeviceSynchronize();
	// exam_over<<<1,1>>>(nparts,cu_overwgt,cu_overthin);
	// cudaDeviceSynchronize();
	printf("val_weight=%d val_thin=%d\n",val_weight,val_thin);

	iter = 0;
	while(val_weight > 0 && val_thin > 0)
	{
		if(iter >= threshold) break;
		iter++;
		// 标记过重边界顶点
		set_bnd_reba<<<(nvtxs + 127) / 128,128>>>(nvtxs,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_where,cu_overwgt,cu_overthin,cu_bnd);

		// 计算顶点移动
		calculate_to_reba<<<nvtxs,32,2 * nparts * sizeof(int)>>>(nvtxs,nparts,graph->cuda_xadj,graph->cuda_adjncy,graph->cuda_adjwgt,\
			graph->cuda_where,cu_bnd,cu_overthin,cu_overwgt,cu_moveto,cu_gain);
		
		// cudaDeviceSynchronize();
		// exam_gain<<<1,1>>>(nvtxs,graph->cuda_where,cu_moveto,cu_gain);
		// cudaDeviceSynchronize();

		init_max<<<(nparts + 31) / 32, 32>>>(nparts,cu_max,cu_maxptr);

		// 选择移动至相同分区的最佳顶点
		calculate_Max_reba<<<(nvtxs + 127) / 128,128>>>(nvtxs,nparts,cu_max,cu_maxptr,cu_moveto,cu_gain);

		// cudaDeviceSynchronize();
		// exam_max<<<1,1>>>(nparts,cu_max,cu_maxptr);
		// cudaDeviceSynchronize();

		// 移动并记录至移动记录数组
		calculate_move<<<(nparts + 31) / 32,32>>>(nparts,cu_maxptr,cu_moveto,graph->cuda_where,cu_moved,graph->cuda_pwgts,graph->cuda_vwgt);

		// calculateSum<<<(nvtxs + 127) / 128,128,nparts * sizeof(int)>>>(nvtxs,nparts,graph->cuda_pwgts,graph->cuda_where,graph->cuda_vwgt);

		calculate_overweight<<<(nparts + 31) / 32, 32>>>(nparts,hunyuangraph_admin->tpwgts[0],graph->tvwgt[0],graph->cuda_pwgts,cu_overwgt,cu_overthin);
		val_weight = thrust::reduce(thrust::device, cu_overwgt, cu_overwgt + nparts);
		val_thin   = thrust::reduce(thrust::device, cu_overthin, cu_overthin + nparts);

		// cudaDeviceSynchronize();
		// exam_pwgts<<<1,1>>>(nparts,hunyuangraph_admin->tpwgts[0],graph->tvwgt[0],graph->cuda_pwgts);
		// cudaDeviceSynchronize();

		// cudaDeviceSynchronize();
		// exam_over<<<1,1>>>(nparts,cu_overwgt,cu_overthin);
		// cudaDeviceSynchronize();

		// printf("val_weight=%d val_thin=%d\n",val_weight,val_thin);
	}

	printf("iter=%d\n",iter);

	cudaFree(cu_bnd);
	cudaFree(cu_moved);
	cudaFree(cu_moveto);
	cudaFree(cu_gain);
	cudaFree(cu_max);
	cudaFree(cu_maxptr);
	cudaFree(cu_overwgt);
	cudaFree(cu_overthin);
}

/*CUDA-kway parjection*/
__global__ void projectback(int *where, int *cwhere, int *cmap, int nvtxs)
{
  int ii = blockIdx.x * blockDim.x + threadIdx.x;

  if(ii < nvtxs)
  {
    int t = cmap[ii];
    where[ii] = cwhere[t];
  }
}

/*Kway parjection*/
void hunyuangraph_kway_project(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph)
{       
  int nvtxs = graph->nvtxs;
  hunyuangraph_graph_t *cgraph = graph->coarser;
  
  projectback<<<(nvtxs + 127) / 128,128>>>(graph->cuda_where,cgraph->cuda_where,graph->cuda_cmap,nvtxs);
}

/*Free graph uncoarsening phase params*/
void hunyuangraph_free_uncoarsen(hunyuangraph_admin_t *hunyuangraph_admin,hunyuangraph_graph_t *graph)
{
	lfree_with_check(sizeof(int) * hunyuangraph_admin->nparts * 2,"cuda_que");			//cuda_que
	lfree_with_check(sizeof(int) * 2,"cuda_csr");										//cuda_csr
	// lfree_with_check(sizeof(int) * graph->nvtxs,"cuda_g");								//cuda_g
	lfree_with_check(sizeof(int) * graph->nvtxs,"cuda_bt");								//cuda_bt
	lfree_with_check(sizeof(int) * graph->nvtxs,"cuda_bn");								//cuda_bn
	lfree_with_check(sizeof(int),"bndnum");												//bndnum
	lfree_with_check(sizeof(int) * hunyuangraph_admin->nparts,"minwgt");				//minwgt
	lfree_with_check(sizeof(int) * hunyuangraph_admin->nparts,"maxwgt");				//maxwgt
	lfree_with_check(sizeof(int) * hunyuangraph_admin->nparts,"tpwgts");				//tpwgts
	lfree_with_check(sizeof(int) * hunyuangraph_admin->nparts,"pwgts");					//pwgts
	lfree_with_check(sizeof(int) * graph->nvtxs,"bnd");									//bnd
	lfree_with_check(sizeof(int) * graph->nvtxs,"where");								//where
	if(graph->cuda_cmap != NULL) lfree_with_check(sizeof(int) * graph->nvtxs,"cmap");	//cmap;
	lfree_with_check(sizeof(int) * graph->nedges,"adjwgt");								//adjwgt
	lfree_with_check(sizeof(int) * graph->nedges,"adjncy");								//adjncy
	lfree_with_check(sizeof(int) * (graph->nvtxs + 1),"xadj");							//xadj
	lfree_with_check(sizeof(int) * graph->nvtxs,"vwgt");								//vwgt
	// cudaFree(graph->cuda_adjwgt);
	// cudaFree(graph->cuda_adjncy);
	// cudaFree(graph->cuda_xadj);
	// cudaFree(graph->cuda_vwgt);
//   cudaFree(graph->cuda_cmap);
//   cudaFree(graph->cuda_maxwgt);
//   cudaFree(graph->cuda_minwgt);
//   cudaFree(graph->cuda_where);
//   cudaFree(graph->cuda_pwgts);
//   cudaFree(graph->cuda_bnd);
//   cudaFree(graph->cuda_bndnum);
//   cudaFree(graph->cuda_real_bnd_num);
//   cudaFree(graph->cuda_real_bnd);
//   cudaFree(graph->cuda_tpwgts);
}

void hunyuangraph_GPU_uncoarsen(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, hunyuangraph_graph_t *cgraph)
{
    cudaDeviceSynchronize();
	gettimeofday(&begin_part_mallocrefine,NULL);
	Mallocinit_refineinfo(hunyuangraph_admin,cgraph);
	cudaDeviceSynchronize();
	gettimeofday(&end_part_mallocrefine,NULL);
	part_mallocrefine += (end_part_mallocrefine.tv_sec - begin_part_mallocrefine.tv_sec) * 1000 + (end_part_mallocrefine.tv_usec - begin_part_mallocrefine.tv_usec) / 1000.0;

	cudaDeviceSynchronize();
	gettimeofday(&begin_part_krefine,NULL);
	hunyuangraph_k_refinement(hunyuangraph_admin,cgraph);
	// hunyuangraph_k_refinement_me(hunyuangraph_admin,cgraph);
	cudaDeviceSynchronize();
	gettimeofday(&end_part_krefine,NULL);
	part_krefine += (end_part_krefine.tv_sec - begin_part_krefine.tv_sec) * 1000 + (end_part_krefine.tv_usec - begin_part_krefine.tv_usec) / 1000.0;

	// int *temp;
	// cudaMalloc((void**)&temp,sizeof(int) * cgraph->nvtxs);
	// init_0<<<(cgraph->nvtxs + 127) / 128,128>>>(cgraph->nvtxs,temp);
	// calculate_edgecut<<<cgraph->nvtxs,32>>>(cgraph->nvtxs,cgraph->cuda_xadj,cgraph->cuda_adjncy,cgraph->cuda_adjwgt,cgraph->cuda_where,temp);
	// int e = thrust::reduce(thrust::device, temp, temp + cgraph->nvtxs);
	// printf("first_krefine me edgecut=%d\n",e / 2);
	// cudaFree(temp);

    for(int i=0;;i++)
	{
    	if(cgraph!=graph)
		{
			cgraph = cgraph->finer;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_mallocrefine,NULL);
			// cudaMalloc((void**)&cgraph->cuda_where, cgraph->nvtxs * sizeof(int));
			cudaDeviceSynchronize();
			gettimeofday(&end_part_mallocrefine,NULL);
			part_mallocrefine += (end_part_mallocrefine.tv_sec - begin_part_mallocrefine.tv_sec) * 1000 + (end_part_mallocrefine.tv_usec - begin_part_mallocrefine.tv_usec) / 1000.0;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_map,NULL);
			hunyuangraph_kway_project(hunyuangraph_admin,cgraph);
			cudaDeviceSynchronize();
			gettimeofday(&end_part_map,NULL);
			part_map += (end_part_map.tv_sec - begin_part_map.tv_sec) * 1000 + (end_part_map.tv_usec - begin_part_map.tv_usec) / 1000.0;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_mallocrefine,NULL);
			hunyuangraph_malloc_refineinfo(hunyuangraph_admin,cgraph);   
			cudaDeviceSynchronize();
			gettimeofday(&end_part_mallocrefine,NULL);
			part_mallocrefine += (end_part_mallocrefine.tv_sec - begin_part_mallocrefine.tv_sec) * 1000 + (end_part_mallocrefine.tv_usec - begin_part_mallocrefine.tv_usec) / 1000.0;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_krefine,NULL);
			hunyuangraph_k_refinement(hunyuangraph_admin,cgraph);
			// hunyuangraph_k_refinement_me(hunyuangraph_admin,cgraph);
			cudaDeviceSynchronize();
			gettimeofday(&end_part_krefine,NULL);
			part_krefine += (end_part_krefine.tv_sec - begin_part_krefine.tv_sec) * 1000 + (end_part_krefine.tv_usec - begin_part_krefine.tv_usec) / 1000.0;

			// cudaMalloc((void**)&temp,sizeof(int) * cgraph->nvtxs);
			// init_0<<<(cgraph->nvtxs + 127) / 128,128>>>(cgraph->nvtxs,temp);
			// calculate_edgecut<<<cgraph->nvtxs,32>>>(cgraph->nvtxs,cgraph->cuda_xadj,cgraph->cuda_adjncy,cgraph->cuda_adjwgt,cgraph->cuda_where,temp);
			// int e = thrust::reduce(thrust::device, temp, temp + cgraph->nvtxs);
			// printf("first_krefine me edgecut=%d\n",e / 2);
			// cudaFree(temp);

			// printf("pwgts:\n");
			// cudaDeviceSynchronize();
			// exam_que<<<1,1>>>(hunyuangraph_admin->nparts-1,cgraph->cuda_pwgts);
			// cudaDeviceSynchronize();

			hunyuangraph_free_uncoarsen(hunyuangraph_admin,cgraph->coarser);
    	} 
    	else 
      		break; 
  	}
}

/*Graph kway-partition algorithm*/
void hunyuangraph_kway_partition(hunyuangraph_admin_t *hunyuangraph_admin, hunyuangraph_graph_t *graph, int *part)
{
	hunyuangraph_graph_t *cgraph;

	cudaDeviceSynchronize();
	gettimeofday(&begin_part_coarsen,NULL);
	cgraph = hunyuangarph_coarsen(hunyuangraph_admin, graph);
	cudaDeviceSynchronize();
	gettimeofday(&end_part_coarsen,NULL);
	part_coarsen += (end_part_coarsen.tv_sec - begin_part_coarsen.tv_sec) * 1000 + (end_part_coarsen.tv_usec - begin_part_coarsen.tv_usec) / 1000.0;

	// printf("Coarsen end:cnvtxs=%d cnedges=%d\n",cgraph->nvtxs,cgraph->nedges);

	cudaDeviceSynchronize();
	gettimeofday(&begin_part_init,NULL);
	hunyuangarph_initialpartition(hunyuangraph_admin, cgraph);
	cudaDeviceSynchronize();
	gettimeofday(&end_part_init,NULL);
	part_init += (end_part_init.tv_sec - begin_part_init.tv_sec) * 1000 + (end_part_init.tv_usec - begin_part_init.tv_usec) / 1000.0;
	
	// hunyuangraph_memcpy_coarsentoinit(cgraph);
	// int e = hunyuangraph_computecut(cgraph, cgraph->where);
	// printf("edgecut=%d\n",e);
	// printf("Init partition end\n");

	cudaDeviceSynchronize();
	gettimeofday(&begin_part_uncoarsen,NULL);
    hunyuangraph_GPU_uncoarsen(hunyuangraph_admin, graph, cgraph);
	/*cudaDeviceSynchronize();
	gettimeofday(&begin_part_mallocrefine,NULL);
	Mallocinit_refineinfo(hunyuangraph_admin,cgraph);
	cudaDeviceSynchronize();
	gettimeofday(&end_part_mallocrefine,NULL);
	part_mallocrefine += (end_part_mallocrefine.tv_sec - begin_part_mallocrefine.tv_sec) * 1000 + (end_part_mallocrefine.tv_usec - begin_part_mallocrefine.tv_usec) / 1000.0;

	cudaDeviceSynchronize();
	gettimeofday(&begin_part_krefine,NULL);
	hunyuangraph_k_refinement(hunyuangraph_admin,cgraph);
	cudaDeviceSynchronize();
	gettimeofday(&end_part_krefine,NULL);
	part_krefine += (end_part_krefine.tv_sec - begin_part_krefine.tv_sec) * 1000 + (end_part_krefine.tv_usec - begin_part_krefine.tv_usec) / 1000.0;
  
  	for(int i=0;;i++)
	{
    	if(cgraph!=graph)
		{
			cgraph = cgraph->finer;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_mallocrefine,NULL);
			cudaMalloc((void**)&cgraph->cuda_where, cgraph->nvtxs * sizeof(int));
			cudaDeviceSynchronize();
			gettimeofday(&end_part_mallocrefine,NULL);
			part_mallocrefine += (end_part_mallocrefine.tv_sec - begin_part_mallocrefine.tv_sec) * 1000 + (end_part_mallocrefine.tv_usec - begin_part_mallocrefine.tv_usec) / 1000.0;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_map,NULL);
			hunyuangraph_kway_project(hunyuangraph_admin,cgraph);
			cudaDeviceSynchronize();
			gettimeofday(&end_part_map,NULL);
			part_map += (end_part_map.tv_sec - begin_part_map.tv_sec) * 1000 + (end_part_map.tv_usec - begin_part_map.tv_usec) / 1000.0;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_mallocrefine,NULL);
			hunyuangraph_malloc_refineinfo(hunyuangraph_admin,cgraph);   
			cudaDeviceSynchronize();
			gettimeofday(&end_part_mallocrefine,NULL);
			part_mallocrefine += (end_part_mallocrefine.tv_sec - begin_part_mallocrefine.tv_sec) * 1000 + (end_part_mallocrefine.tv_usec - begin_part_mallocrefine.tv_usec) / 1000.0;

			cudaDeviceSynchronize();
			gettimeofday(&begin_part_krefine,NULL);
			hunyuangraph_k_refinement(hunyuangraph_admin,cgraph);
			cudaDeviceSynchronize();
			gettimeofday(&end_part_krefine,NULL);
			part_krefine += (end_part_krefine.tv_sec - begin_part_krefine.tv_sec) * 1000 + (end_part_krefine.tv_usec - begin_part_krefine.tv_usec) / 1000.0;

			hunyuangraph_free_uncoarsen(cgraph->coarser);
    	} 
    	else 
      		break; 
  	}*/
	cudaDeviceSynchronize();
	gettimeofday(&end_part_uncoarsen,NULL);
	part_uncoarsen += (end_part_uncoarsen.tv_sec - begin_part_uncoarsen.tv_sec) * 1000 + (end_part_uncoarsen.tv_usec - begin_part_uncoarsen.tv_usec) / 1000.0;

	// printf("Uncoarsen end\n");
}

/*Graph partition algorithm*/
void hunyuangraph_PartGraph(int *nvtxs,  int *xadj, int *adjncy, int *vwgt,int *adjwgt, \
  int *nparts, float *tpwgts, float *ubvec, int *part)
{
    hunyuangraph_graph_t *graph;
    hunyuangraph_admin_t *hunyuangraph_admin;

    hunyuangraph_admin = hunyuangraph_set_graph_admin(*nparts, tpwgts, ubvec);

    graph = hunyuangraph_set_first_level_graph(*nvtxs, xadj, adjncy, vwgt, adjwgt);

    hunyuangraph_set_kway_bal(hunyuangraph_admin, graph);

    hunyuangraph_admin->Coarsen_threshold = hunyuangraph_max((*nvtxs) / (20 * (hunyuangraph_compute_log2(*nparts))),30 * (*nparts));
    
    hunyuangraph_admin->nIparts = (hunyuangraph_admin->Coarsen_threshold == 30 * (*nparts) ? 4 : 5);

	Malloc_GPU_Memory(graph->nvtxs,graph->nedges);

	// cudaMalloc((void**)&cu_bn, graph->nvtxs * sizeof(int));
	// cudaMalloc((void**)&cu_bt, graph->nvtxs * sizeof(int));
	// cudaMalloc((void**)&cu_g,  graph->nvtxs * sizeof(int));
	// cudaMalloc((void**)&cu_csr, 2 * sizeof(int));
	// cudaMalloc((void**)&cu_que, 2 * hunyuangraph_admin->nparts * sizeof(int));

	// cu_bn  = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_bn");
	// cu_bt  = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_bt");
	// cu_g   = (int *)lmalloc_with_check(sizeof(int) * graph->nvtxs,"cu_g");
	// cu_csr = (int *)lmalloc_with_check(sizeof(int) * 2,"cu_csr");
	// cu_que = (int *)lmalloc_with_check(sizeof(int) * hunyuangraph_admin->nparts * 2,"cu_que");

	hunyuangraph_malloc_original_coarseninfo(hunyuangraph_admin,graph);

	// printf("begin partition\n");
	cudaDeviceSynchronize();
	gettimeofday(&begin_part_all,NULL);
	hunyuangraph_kway_partition(hunyuangraph_admin,graph,part);
	cudaDeviceSynchronize();
	gettimeofday(&end_part_all,NULL);
	part_all += (end_part_all.tv_sec - begin_part_all.tv_sec) * 1000 + (end_part_all.tv_usec - begin_part_all.tv_usec) / 1000.0;
	// printf("end partition\n");

	cudaMemcpy(part,graph->cuda_where, graph->nvtxs * sizeof(int), cudaMemcpyDeviceToHost);

	hunyuangraph_free_uncoarsen(hunyuangraph_admin,graph);

	// lfree_with_check(sizeof(int) * hunyuangraph_admin->nparts * 2,"cu_que");	//cu_que
	// lfree_with_check(sizeof(int) * 2,"cu_csr");									//cu_csr
	// lfree_with_check(sizeof(int) * graph->nvtxs,"cu_g");						//cu_g
	// lfree_with_check(sizeof(int) * graph->nvtxs,"cu_bt");						//cu_bt
	// lfree_with_check(sizeof(int) * graph->nvtxs,"cu_bn");						//cu_bn

	Free_GPU_Memory();

	printf("Max memory used of CPU: %10ldB %10ldKB %10ldMB %10ldGB\n",used_by_init_max,used_by_init_max / 1024,\
		used_by_init_max / 1024 / 1024, used_by_init_max/1024 / 1024 / 1024);
}

/*Error exit*/
void hunyuangraph_error_exit(char *f_str,...)
{
  va_list a;
  va_start(a,f_str);
  vfprintf(stderr,f_str,a);
  va_end(a);

  if (strlen(f_str)==0||f_str[strlen(f_str)-1]!='\n'){
    fprintf(stderr,"\n");
  }

  fflush(stderr);

  if(1)
    exit(-2);
}

/*Open file*/
FILE *hunyuangraph_fopen(char *fname, char *mode, const char *msg)
{
  FILE *fp;
  char error_message[8192];
  fp=fopen(fname, mode);
  if(fp!=NULL){
    return fp;
  }
  sprintf(error_message,"file: %s, mode: %s, [%s]",fname,mode,msg);
  perror(error_message);
  hunyuangraph_error_exit("Failed on file fopen()\n");
  return NULL;
}

/*Read graph file*/
hunyuangraph_graph_t *hunyuangraph_readgraph(char *filename)
{
  int i,k,fmt,nfields,readew,readvw,readvs,edge,ewgt;
  int *xadj,*adjncy,*vwgt,*adjwgt;
  char *line=NULL,fmtstr[256],*curstr,*newstr;
  size_t lnlen=0;
  FILE *fpin;

  hunyuangraph_graph_t *graph;
  graph = hunyuangraph_create_cpu_graph();

  fpin = hunyuangraph_fopen(filename,"r","Readgraph: Graph");

  do{
    if(getline(&line,&lnlen,fpin)==-1){ 
      hunyuangraph_error_exit("Premature end of input file: file: %s\n", filename);
    }
  }while(line[0]=='%');

  fmt= 0;
  nfields = sscanf(line, "%d %d %d", &(graph->nvtxs), &(graph->nedges), &fmt);

  if(nfields<2){
    hunyuangraph_error_exit("The input file does not specify the number of vertices and edges.\n");
  }

  if(graph->nvtxs<=0||graph->nedges<=0){
   hunyuangraph_error_exit("The supplied nvtxs:%d and nedges:%d must be positive.\n",graph->nvtxs,graph->nedges);
  }

  if(fmt>111){ 
    hunyuangraph_error_exit("Cannot read this type of file format [fmt=%d]!\n",fmt);
  }

  sprintf(fmtstr,"%03d",fmt%1000);
  readvs=(fmtstr[0]=='1');
  readvw=(fmtstr[1]=='1');
  readew=(fmtstr[2]=='1');

  graph->nedges *=2;

  xadj=graph->xadj=(int*)malloc(sizeof(int)*(graph->nvtxs+1));
  for(i=0;i<graph->nvtxs+1;i++){
    xadj[i]=graph->xadj[i]=0;
  }

  adjncy=graph->adjncy=(int*)malloc(sizeof(int)*(graph->nedges));

  vwgt=graph->vwgt= (int*)malloc(sizeof(int)*(graph->nvtxs));

  for(i=0;i<graph->nvtxs;i++){
    vwgt[i]=graph->vwgt[i]=1;
  }

  adjwgt = graph->adjwgt=(int*)malloc(sizeof(int)*(graph->nedges));
  for(i=0;i<graph->nedges;i++){
    adjwgt[i]=graph->adjwgt[i]=1;
  }

  for(xadj[0]=0,k=0,i=0;i<graph->nvtxs;i++){
    do{
      if(getline(&line,&lnlen,fpin)==-1){
      hunyuangraph_error_exit("Premature end of input file while reading vertex %d.\n", i+1);
      } 
    }while(line[0]=='%');

    curstr=line;
    newstr=NULL;

    if(readvw){
      vwgt[i]=strtol(curstr, &newstr, 10);

      if(newstr==curstr){
        hunyuangraph_error_exit("The line for vertex %d does not have enough weights "
          "for the %d constraints.\n", i+1, 1);
      }
      if(vwgt[i]<0){
        hunyuangraph_error_exit("The weight vertex %d and constraint %d must be >= 0\n", i+1, 0);
      }
      curstr = newstr;
    }

    while(1){
      edge=strtol(curstr,&newstr,10);
      if(newstr==curstr){
        break; 
      }

      curstr=newstr;
      if (edge< 1||edge>graph->nvtxs){
        hunyuangraph_error_exit("Edge %d for vertex %d is out of bounds\n",edge,i+1);
      }

      ewgt=1;

      if(readew){
        ewgt=strtol(curstr,&newstr,10);

        if(newstr==curstr){
          hunyuangraph_error_exit("Premature end of line for vertex %d\n", i+1);
        }

        if(ewgt<=0){
          hunyuangraph_error_exit("The weight (%d) for edge (%d, %d) must be positive.\n",    ewgt, i+1, edge);
        }

        curstr=newstr;
      }

      if(k==graph->nedges){
        hunyuangraph_error_exit("There are more edges in the file than the %d specified.\n", graph->nedges/2);
      }

      adjncy[k]=edge-1;
      adjwgt[k]=ewgt;
      k++;

    } 
    xadj[i+1]=k;

  }
  fclose(fpin);

  if(k!=graph->nedges){
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("In the first line of the file, you specified that the graph contained\n"
      "%d edges. However, I only found %d edges in the file.\n", graph->nedges/2,k/2);
    if(2*k==graph->nedges){
      printf("\n *> I detected that you specified twice the number of edges that you have in\n");
      printf("    the file. Remember that the number of edges specified in the first line\n");
      printf("    counts each edge between vertices v and u only once.\n\n");
    }
    printf("Please specify the correct number of edges in the first line of the file.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }
  free(line);
  return graph;
}

/*Write to file*/
void hunyuangraph_writetofile(char *fname, int *part, int n, int nparts)
{
  FILE *fpout;
  int i;
  char filename[1280000];
  sprintf(filename, "%s.part.%d", fname, nparts);

  fpout = hunyuangraph_fopen(filename, "w", __func__);

  for (i=0; i<n; i++){
    fprintf(fpout,"%d\n",part[i]);
  }

  fclose(fpout);
}

/*Main function*/
int main(int argc, char **argv)
{  
    cudaSetDevice(0);

    char *filename = (argv[1]);
    int nparts     = atoi(argv[2]);

    hunyuangraph_graph_t *graph = hunyuangraph_readgraph(filename);

    printf("graph:%s %d %d\n",filename,graph->nvtxs,graph->nedges);
  
    int *part = (int*)malloc(sizeof(int) * graph->nvtxs);

    float tpwgts[nparts];
    for(int i = 0;i < nparts;i++)
        tpwgts[i] = 1.0 / nparts;
  
    float ubvec=1.03;
    
    hunyuangraph_PartGraph(&graph->nvtxs, graph->xadj, graph->adjncy, graph->vwgt, graph->adjwgt, &nparts, tpwgts, &ubvec, part);

    int e = hunyuangraph_computecut(graph, part);

	printf("Hunyuangraph-Partition-end\n");
	printf("Hunyuangraph_Partition_time= %10.2lf ms\n",part_all);
	printf("------Coarsen_time=          %10.2lf ms\n",part_coarsen);
	printf("------Init_time=             %10.2lf ms\n",part_init);
	printf("------Uncoarsen_time=        %10.2lf ms\n",part_uncoarsen);
	printf("------else_time=             %10.2lf ms\n",part_all - (part_coarsen + part_init + part_uncoarsen));
	printf("edge-cut=                    %10d\n",e);

	printf("\n");
    printf("\n");
}