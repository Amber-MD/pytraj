html:
	python ./scripts/check_libs.py
	(cd doc && make html)
	python ./scripts/assert_no_errors.py

clean:
	git clean -fxd

push:
	sh ./doc/scripts/add_git_and_push.sh

pull:
	git pull upstream gh-pages
