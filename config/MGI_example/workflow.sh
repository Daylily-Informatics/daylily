#!/bin/bash
 
set -e
err_report() {
  echo "Error on line $1, error code:$?"
}
trap 'err_report $LINENO' ERR
 
 
workflow_dir=/path/to/test-project
cd $workflow_dir
 
# if you need to use a mysql database
# mysql_dir=$workflow_dir/mysql-directory
# mkdir -p $mysql_dir/run/mysqld
# mkdir -p $mysql_dir/lib/mysql
# mkdir -p $mysql_dir/log/mysql
# mysql_install_db --user=$USER --basedir=/usr/ --ldata=$mysql_dir/lib/mysql/
# mysqld_safe --defaults-file=/$mysql_dir/my.cnf &
# for i in `seq 1 28`  ; do 
#     if [ ! -e "/tmp/mysqld.sock" ] ; then 
#         sleep $i; 
#     fi; 
# done

# echo "create database cromwell; create user 'cromwell'@'localhost' identified by 'test4cromwell'; grant all privileges on *.* to 'cromwell'@localhost;" \
# | mysql -u root --socket=/tmp/mysqld.sock
 
java -Dconfig.file=config/MGI_example/app.conf -jar /usr/local/bin/cromwell.jar run config/MGI_example/example.wdl
#/usr/bin/java -Dconfig.file=$workflow_dir/application.conf -jar /cromwell/cromwell.jar run $workflow_dir/example.wdl -
 
# If you setup a mysql database, make sure you properly shut it down
# /usr/bin/mysqladmin -u root --socket /tmp/mysqld.sock shutdown