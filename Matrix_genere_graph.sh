input="graph.csv"

{
  read
  i=1
  while IFS=',' read -r  Name 
  do
    ./Matrix_genere_graph mtx/$Name.mtx > graph/$Name.graph
 
    i=`expr $i + 1`
  done 
} < "$input"