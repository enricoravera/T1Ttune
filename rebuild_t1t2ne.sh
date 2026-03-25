#! /bin/bash

###----------------------------------------------------------------------------------------------------------------------------------------------------------
### Use option -d to rebuild the documentation as well (both html and latex).
###
###----------------------------------------------------------------------------------------------------------------------------------------------------------

# ! # EDIT HERE # ! #
# -------------------
source_dir="/data_old/NMR/HetRelax/t1tune_main"     # Location of the package source code
version="0.0.2.dev1"
# -------------------

redir_file="rebuild_t1t2ne"		# Name of the redirection files .log and .err
mydir=$(pwd)						# Where am I right now


# Check existence of log file: if there already is, remove it to allow fresh writing
if [ -f redir_file.log ]; then
	rm redir_file.log
	rm redir_file.err
fi

cd $source_dir		# Go to the source directory

# Replace the versions
#	LICENSE
sed -i "s/Version: \"[^\"]*\"/Version = $version/g" LICENSE.txt
#	init.py
sed -i "s/__version__ = \"[^\"]*\"/__version__ = \"$version\"/g" ./t1t2ne/__init__.py
#	conf.py
sed -i "s/release = \"[^\"]*\"/release = \"$version\"/g" ./docs/source/conf.py
#	index.rst
sed -i "s/\*\*Version\*\*: .*/**Version**: $version/g" ./docs/source/index.rst
#       pyproject.toml
sed -i "s/version = \"[^\"]*\"/version = \"$version\"/g" ./pyproject.toml

# Rebuild the package
echo "Rebuilding t1t2ne version $version..."
python -m build ${source_dir} >> ${mydir}/${redir_file}.log 2>> ${mydir}/${redir_file}.err

echo 'Done.'

# Install the package
echo 'Reinstalling t1t2ne...'
pip install -e ${source_dir} >> ${mydir}/${redir_file}.log 2>> ${mydir}/${redir_file}.err
echo 'Done.'


while getopts "d" opt; do
	case $opt in
		d)
			# Rebuild documentation
			echo 'Rebuilding documentation...'
			cd ${source_dir}/docs
			make clean >> ${mydir}/${redir_file}.log 2>> ${mydir}/${redir_file}.err
			make html >> ${mydir}/${redir_file}.log 2>> ${mydir}/${redir_file}.err
			make latexpdf >> ${mydir}/${redir_file}.log 2>> ${mydir}/${redir_file}.err
			echo 'Done.'
			;;
		\?)
			echo "Opzione non valida: -$OPTARG" >&2
			exit 1
			;;
	esac
done


# Go back to the original directory
cd $mydir

# If I'm here it means all ran smoothly. Remove the .log and .err files.
rm ${mydir}/${redir_file}.log
rm ${mydir}/${redir_file}.err


