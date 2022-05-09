#devtools::create('~/nfs/rcode/visutils/')
devtools::document()
devtools::load_all()
?log10p1

git add -u
git commit -m "updated documentation"
git push -u origin main
