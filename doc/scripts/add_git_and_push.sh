#!/bin/sh

git ls-files --modified | xargs git add
git commit -m 'update doc'
git push upstream gh-pages
