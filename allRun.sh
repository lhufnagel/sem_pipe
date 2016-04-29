rm -f pipe.sch
rm -f SESSION.NAME
rm -f pipe0.f*

echo pipe    >  SESSION.NAME
echo `readlink -f \`pwd\``/ >>  SESSION.NAME
#
mpiexec -np $1 ./nek5000
#
