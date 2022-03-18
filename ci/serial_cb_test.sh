#!/bin/bash


SELF_BUILD_DIR=/build
SELF_INSTALL_DIR=/opt/self

moveGCNOFiles(){

  for file in $(ls $SELF_BUILD_DIR/src/*.gcno); do
    dest=$(echo $file | sed 's/\.f//g')
    cp $file $dest
  done

  for file in $(ls $SELF_BUILD_DIR/examples/*.gcno); do
    dest=$(echo $file | sed 's/\.f90//g')
    cp $file $dest
  done

}

moveGCDAFiles(){
  # This method is used to apply the correct file extension
  # for gcovr to be able to match up the gcda file to the 
  # source code.
  #
  # This is only needed because the SELF build system applies
  # the .f.o extension to object files created from our .f90
  # files.

  for file in $(ls $SELF_BUILD_DIR/src/*.gcda); do
    dest=$(echo $file | sed 's/\.f//g')
    cp $file $dest
  done

  for file in $(ls $SELF_BUILD_DIR/examples/*.gcda); do
    dest=$(echo $file | sed 's/\.f90//g')
    cp $file $dest
  done

}

cleanupGCDAFiles(){
  rm $SELF_BUILD_DIR/examples/*.gcda
  rm $SELF_BUILD_DIR/src/*.gcda
}

moveGCNOFiles
tomerge=""
for file in $(ls $SELF_INSTALL_DIR/examples/*); do

  echo "Running Test : $file"
  $file
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "$file failed!"
    exit $rc
  fi

  moveGCDAFiles

  # Create a coverage file for this test.
  covFile=$(echo $file | awk -F "/" '{print $NF}')
  gcovr -r /build --json -o "/build/$covFile.json"
  tomerge+="-a /build/$covFile.json "

  cleanupGCDAFiles

done

  gcovr $tomerge
  gcovr $tomerge -x -o /build/coverage.xml
