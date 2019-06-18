
sed "s/!TAG!/${BUILD_TAG}/g" tools/scripts/batch-daint-${partition}.sh.tmpl > tools/scripts/batch_${BUILD_TAG}.sh

sbatch --wait tools/scripts/batch_${BUILD_TAG}.sh
ret=$?

echo ""

echo "----- Result: -----"
if [ "$ret" -eq "0" ]; then
  echo "Test passed."
else
  echo "Test FAILED. (Error: $ret)"
fi

echo ""

ARCHIVE=/project/csstaff/rasolca/jenkins/output/${J_PROJ}/${BUILD_ID}/${BUILD_TAG}
mkdir -p -- $ARCHIVE

echo "----- Output: -----"
cp ${BUILD_TAG}.out ${ARCHIVE}/out
cat ${BUILD_TAG}.out

echo ""

echo "----- Error:  -----"
cp ${BUILD_TAG}.err ${ARCHIVE}/err
cat ${BUILD_TAG}.err

echo ""

echo "----- Log:    -----"
cp build/Testing/Temporary/LastTest.log* ${ARCHIVE}/log
cat build/Testing/Temporary/LastTest.log*

rm -rf $WORKSPACE/build

if [ "$ret" -ne "0" ]; then
  exit 2
fi
