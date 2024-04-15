all:
	gcc Matrix_genere_graph.c -o Matrix_genere_graph
	nvcc -std=c++11 -gencode arch=compute_89,code=sm_89 -O3  Hunyuan.cu -o  Hunyuan  --expt-relaxed-constexpr -w