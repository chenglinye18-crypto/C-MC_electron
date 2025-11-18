/* --------------------- Macros for Range Handling --------------------- */
// These macros ensure indices are within the local range,
//  including ghost cell regions. If an index `i` is outside 
//  the specified range, it's clamped to the nearest boundary.
/**
 * 
 */
#define P_LOCAL_LRANGE(i) (((i) < p_jbegin) ? p_jbegin : i)
#define P_LOCAL_RRANGE(i) (((i) > p_jend) ? p_jend : i)

/**
 * @brief
 * Ensures that the index `i` is not less than `p_jbegin_ghost`
 *  If i is less than p_jbegin_ghost, it returns p_jbegin_ghost; 
 *  otherwise, it returns i.
 */
#define P_LOCAL_LRANGE_GHOST(i) (((i) < p_jbegin_ghost) ? p_jbegin_ghost : i)

/**
 * @brief
 * Ensures that the index i is not greater than p_jend_ghost. 
 *  If i is greater than p_jend_ghost, it returns p_jend_ghost; 
 *  otherwise, it returns i.
 */
#define P_LOCAL_RRANGE_GHOST(i) (((i) > p_jend_ghost) ? p_jend_ghost : i)


/* --------------------- Macros for Index Calculation --------------------- */
// These macros calculate the 2D index for the given coordinate `(i,k)` 
//  in a 3D grid.
// 下面两个宏看起来分别是计算 Z 方向的某个cell 或 point
/**
 * @brief 统计 (i,k) 点处 xz 面上的 cell 数 和 Point 数
 * 
 * @details
 * c_numz 代表了 Device mesh 中，Z方向上一行的 cell 数
 * i 为 x 方向上的坐标
 * 因此， i*c_numz 为到 点(i,k)处，沿着z轴，i行的 cell 数目
 * 又因为，k 为 z 方向上的坐标，
 *  因此 i*c_numz + k 代表了 到点(i,k)处，xz 面上的 cell 总数
 * 
 * 想不明白就画个图
 */
#define C_2DINDEX(i,k) ((i) * c_numz + (k))
#define P_2DINDEX(i,k) ((i) * p_numz + (k))

/* --------------------- Macros for Linear Indexing --------------------- */
// These macros convert 3D indices `(i,j,k)` into linear indices for accessing
//  arrays, considering both regular and ghost cells.
// "L" means Local while "G" means Global.
/**
 * @brief
 * Converts 3D indices (i, j, k) and additional indices `idx` and `tot` into a single linear index 
 *      for a 3D array with ghost cells. 
 *  `c_numxz` is the product of the dimensions in the x and z directions.
 * 
 * 
 */
#define C_LINDEX_GHOST(i,j,k,idx,tot) ((((j) - c_jbegin_ghost) * c_numxz + C_2DINDEX(i,k)) * (tot) + (idx))
#define C_LINDEX(i,j,k,idx,tot) ((((j) - c_jbegin) * c_numxz + C_2DINDEX(i,k)) * (tot) + (idx))
#define P_LINDEX(i,j,k,idx,tot) ((((j) - p_jbegin) * p_numxz + P_2DINDEX(i,k)) * (tot) + (idx))
#define P_LINDEX_GHOST(i,j,k,idx,tot) ((((j) - p_jbegin_ghost) * p_numxz + P_2DINDEX(i,k)) * (tot) + (idx))

/**
 *  @brief 
 * 取得位置为 (i,j,k) 的 ghost cell 对应的 local index
 */
#define C_LINDEX_GHOST_ONE(i,j,k) C_LINDEX_GHOST(i,j,k,0,1)

/**
 * @brief 
 *  到点 (i,j,k) 处，cell 的个数，也可以理解成 索引值
 * @details
 *  等价于
 *  ( (j) - c_jbegin) * c_numxz + C_2DINDEX(i,k) ) * 1 + 0
 *  其中，j-c_jbegin 代表了到第j个cell处，y方向上 cell 个数
 *      因此，(j) - c_jbegin) * c_numxz 代表了到第 j 个 cell 处，cell 的总和    
 *  而 C_2DINDEX(i,k) 代表此时在 xz 方向上的 cell 数目
 * 
 *  因此，C_LINDEX_ONE(i,j,k) 代表了 到点 (i,j,k) 处，cell 的个数
 *  
 * @attention
 * c_jbegin 代表了当前 Processor 中沿 y 方向的 起始cell 的global index
 *  而不是整个 device 的
 * 因此这个 macro 是用于计算 local index
 */
#define C_LINDEX_ONE(i,j,k) C_LINDEX(i,j,k,0,1)
#define P_LINDEX_ONE(i,j,k) P_LINDEX(i,j,k,0,1)
#define P_LINDEX_ONE_GHOST(i,j,k) P_LINDEX_GHOST(i,j,k,0,1)
#define P_LINDEX_N_GHOST(i,j,k,jb) (((j) - jb) * p_numxz + P_2DINDEX(i,k))

/**
 * @brief
 * 到点(i,j,k)处，cell 或 Point 的总数
 * 或者说，计算的结果为一个 global index (但要小心 ±1 的索引问题)
 * 看起来返回的应该是 global index
 * 
 * @details 
 *  P_2DINDEX 用于计算点(i,k)上
 *  j 为 y 方向上的坐标
 *  (j) * c_numxz 代表了前面 j 行上的 cell 总数
 *  那么，整个 Macro 代表了 到点(i,j,k)处，cell 的总数
 * 
 * 想不明白就画个图
 */
#define C_GINDEX(i,j,k) ((j) * c_numxz + C_2DINDEX(i,k))
#define P_GINDEX(i,j,k) ((j) * p_numxz + P_2DINDEX(i,k))

#define P_SHIFT_NONOVERLAP_Y(j) ((j) - p_jbegin_nonoverlap)

#define PAR_EXIST 0x1
#define PAR_UNFINISHED 0x2

#define PAR_DBL_NUM 9
#define PAR_INT_NUM 9

#define P_QCLINDEX(i,j,k) (((j) - p_qc_jbegin) * p_numxz + P_2DINDEX(i,k))

#define DEBUG_DATA_PRINT
#define DEBUG_PRINT 100

#define NOT_GHOST_BC(j) (((j) != 0) && ((j) != p_numy - 1)) 
#define NOT_GHOST_CELL(j) (((j) != 0) && ((j) != c_numy - 1)) 

#define MY_ZERO 1e-13

#define BETWEEN01(x) (((x) <= 1.0 + MY_ZERO) && ((x) >= - MY_ZERO))

#define NODE_SILICON (1)
#define NODE_OXIDE (2)
#define NODE_UP_BC (1 << 2)
#define NODE_DOWN_BC (1 << 3)
#define NODE_FRONT_BC (1 << 4)
#define NODE_BACK_BC (1 << 5)
#define NODE_LEFT_BC (1 << 6)
#define NODE_RIGHT_BC (1 << 7)
#define NODE_QUANTUM (1 << 8)
#define STARS "**********"
