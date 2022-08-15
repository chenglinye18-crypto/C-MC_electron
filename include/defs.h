#define P_LOCAL_LRANGE(i) (((i) < p_jbegin) ? p_jbegin : i)
#define P_LOCAL_RRANGE(i) (((i) > p_jend) ? p_jend : i)

#define P_LOCAL_LRANGE_GHOST(i) (((i) < p_jbegin_ghost) ? p_jbegin_ghost : i)
#define P_LOCAL_RRANGE_GHOST(i) (((i) > p_jend_ghost) ? p_jend_ghost : i)

#define C_2DINDEX(i,k) ((i) * c_numz + (k))
#define P_2DINDEX(i,k) ((i) * p_numz + (k))

#define C_LINDEX_GHOST(i,j,k,idx,tot) ((((j) - c_jbegin_ghost) * c_numxz + C_2DINDEX(i,k)) * (tot) + (idx))
#define C_LINDEX(i,j,k,idx,tot) ((((j) - c_jbegin) * c_numxz + C_2DINDEX(i,k)) * (tot) + (idx))
#define P_LINDEX(i,j,k,idx,tot) ((((j) - p_jbegin) * p_numxz + P_2DINDEX(i,k)) * (tot) + (idx))
#define P_LINDEX_GHOST(i,j,k,idx,tot) ((((j) - p_jbegin_ghost) * p_numxz + P_2DINDEX(i,k)) * (tot) + (idx))

// C: Cell, P: Point
#define C_LINDEX_GHOST_ONE(i,j,k) C_LINDEX_GHOST(i,j,k,0,1)
#define C_LINDEX_ONE(i,j,k) C_LINDEX(i,j,k,0,1)
#define P_LINDEX_ONE(i,j,k) P_LINDEX(i,j,k,0,1)
#define P_LINDEX_ONE_GHOST(i,j,k) P_LINDEX_GHOST(i,j,k,0,1)
#define P_LINDEX_N_GHOST(i,j,k,jb) (((j) - jb) * p_numxz + P_2DINDEX(i,k))

#define P_GINDEX(i,j,k) ((j) * p_numxz + P_2DINDEX(i,k))
#define C_GINDEX(i,j,k) ((j) * c_numxz + C_2DINDEX(i,k))

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
