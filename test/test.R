#devtools::create('~/nfs/rcode/visutils/')
devtools::document()
devtools::load_all()
?log10p1

# usethis::create_github_token()
# gitcreds::gitcreds_set()
# usethis::edit_r_environ()
options(timeout=400)
devtools::install_github("iaaka/visutils")

devtools::install_local()
# git add -u
# git commit -m "add spot merging"
# git push -u origin main
?dotPlot


