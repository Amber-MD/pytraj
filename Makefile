html:
	(cd doc && make html)

clean:
	git clean -fxd

push:
	sh ./doc/scripts/add_git_and_push.sh
