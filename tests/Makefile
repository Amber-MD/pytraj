test:
	py.test -v . -n4 --cov=pytraj --cov-report=html

autoflake:
	autoflake -i --remove-all-unused-imports --remove-unused-variables *.py */*.py

clean:
	git clean -fdx .
