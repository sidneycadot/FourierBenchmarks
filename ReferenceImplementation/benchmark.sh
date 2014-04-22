#! /bin/sh

for TIME_LIMIT in 1 2 5 10 20 50 100 200 500 1000
do
    ./Test_FastFourierTransform_ReferenceImplementation $TIME_LIMIT | tee `hostname`_${TIME_LIMIT}ms.txt
done
