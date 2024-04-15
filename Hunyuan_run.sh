input="graph.csv"

{
  read
  i=1
  while IFS=',' read -r  Name 
  do
    ./Hunyuan graph/$Name.graph 8 >> result/$Name\_part8.txt
    ./Hunyuan graph/$Name.graph 32 >> result/$Name\_part32.txt
    ./Hunyuan graph/$Name.graph 64 >> result/$Name\_part64.txt
    ./Hunyuan graph/$Name.graph 128 >> result/$Name\_part128.txt
    ./Hunyuan graph/$Name.graph 256 >> result/$Name\_part256.txt
 
    i=`expr $i + 1`
  done 
} < "$input"