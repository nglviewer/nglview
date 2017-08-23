html:
	cp -rf ../nglview/docs .
	python  ./scripts/fix_readme_link.py ../nglview/README.md
	cp docs/ngl*.{js,html} dev/
	cp docs/ngl*.{js,html} latest
	(cd docs && make html)

clean:
	git clean -fxd

push:
	git push upstream gh-pages

pull:
	git pull upstream gh-pages
