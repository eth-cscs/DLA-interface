
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

echo "----- Output: -----"
cat ${BUILD_TAG}.out

echo ""

echo "----- Error:  -----"
cat ${BUILD_TAG}.err

echo ""

echo "----- Log:    -----"
cp build_${BUILD_TAG}/Testing/Temporary/LastTest.log* /project/csstaff/rasolca/jenkins/output/${BUILD_TAG}-out.log
cat build_${BUILD_TAG}/Testing/Temporary/LastTest.log*

if [ "$ret" -ne "0" ]; then
  exit 2
fi
