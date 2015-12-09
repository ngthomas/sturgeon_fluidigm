
if [ $# -ne 1 ]; then
echo Syntax:
echo $(basename $0) NumProc
echo "
NumProc is the number of processors to spread this job over.
You should set it to be the number of cores on your machine 
or one less than that, etc.

"
exit 1;
fi;

cat ../input/Commands.txt  | ../script/parallel -P$1


