#! /bin/bash

src=$1
dst=$2

# Copy tree $src to tree $dst

# Both $src and $dst must be directories.  $dst is created if it does
# not exist

# All files in the source tree are checked; if they already exist in
# the destination tree and are identical, they are ignored, otherwise
# they are copied.  Missing directories are created.

# All files in the destination tree are checked; if they do not exist
# in the source tree, they are deleted.

if test -z "$src" || test -z "$dst" || test "$src" = "$dst"; then
    echo "Usage: $0 <src> <dst>"
    exit 1
fi
test -d $src || exit 2
test -e $dst || mkdir -p $dst
test -d $dst || exit 3

# Create all directories
for dir in $(cd $src && find . -type d); do
    dstdir=$dst/$dir
    if test -d $dstdir; then
        : # directory exists; do nothing
    else
        echo mkdir $dstdir
        mkdir -p $dstdir
    fi
done

# Delete directories which do not exist
for dir in $(cd $dst && find . -type d); do
    srcdir=$src/$dir
    dstdir=$dst/$dir
    if test -d $srcdir; then
        : # directory exists; do nothing
    else
        echo rm -rf $dstdir
        #rm -rf $dstdir
    fi
done

# Copy files that differ
for file in $(cd $src && find . -type f); do
    srcfile=$src/$file
    dstfile=$dst/$file
    if cmp -s $srcfile $dstfile; then
        : # unchanged; do nothing
    else
        echo cp $srcfile
        cp $srcfile $dstfile
    fi
done

# Delete files which do not exist
for file in $(cd $dst && find . -type f); do
    srcfile=$src/$file
    dstfile=$dst/$file
    if test -e $srcfile; then
        : # file exists; do nothing
    else
        echo rm $dstfile
        rm $dstfile
    fi
done
