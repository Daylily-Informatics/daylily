ls -d1 .tests/unit/*/data/ | parallel 'p=$PWD; cd {} && ln -s ../../../../config config && ln -s ../../../../resources resources ' #&& ln -s ../../../../.test_data .test_data'
cp bin/helpers/common.py .tests/unit/
