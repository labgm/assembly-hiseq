#!/usr/bin/env bash

if [ -f $1 ] ; then
case $1 in
*.tar.bz2) tar -xjvf $1 -O > $2 ;;
*.tar.gz) tar -xzvf $1 -O > $2 ;;
*.bz2) bzip2 -dck $1 > $2 ;;
*.rar) rar p $1 > $2 ;;
*.gz) gunzip -kc $1 > $2 ;;
*.tar) tar -xf $1 -O > $2 ;;
*.tbz2) tar -xjvf $1 -O > $2 ;;
*.tgz) tar -xzvf $1 -O > $2 ;;
*.zip) unzip -p $1 > $2 ;;
*.Z) uncompress -dck $1 > $2 ;;
*.7z) 7z x -so $1 > $2 ;;
*) echo "'$1' could not be extracted" ;;
esac
else
echo "'$1' is not a valid file"
fi
