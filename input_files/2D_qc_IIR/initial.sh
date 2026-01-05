nmol=64
method="did"
xyzpath="/path/to"
mkdir "$nmol"_"$method"

for i in {1..10}
do
mkdir init_water_xyz_nvt
cp "$xyzpath"/"$nmol"_nvt/h2o_"$i"/nvt_1/simu.pos_0.xyz ./
python get_xyz.py
cp -r run_sample_"$nmol"_"$method" "$nmol"_"$method"/"$i"_7e-4
cp -r init_water_xyz_nvt "$nmol"_"$method"/"$i"_7e-4/data
rm -rf init_water_xyz_nvt
rm -rf simu.pos_0.xyz
sed -i "s/newE0=0/newE0=7/" "$nmol"_"$method"/"$i"_7e-4/summit.sh
done
