make MolDyn_NVE
./clean.sh
echo "Run 1"
./MolDyn_NVE
for ((i = 1; i<(($1)); i++))
do
  echo "Run $((i+1))"
  ./rerun.sh
done
python plot.py
python plot_output.py
