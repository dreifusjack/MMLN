### Setting up MMLN Package
## 1. Install Required Packages
install.packages(c("devtools",
                   "usethis",
                   "roxygen2",
                   "testthat"))
#devtools: High-level functions for package development (e.g., devtools::check(), devtools::build(), etc.).
#usethis: Functions to initialize packages, manage R build infrastructure, etc.
#roxygen2: For inline documentation that will generate .Rd files automatically.

## 2. Create a New Package Structure
usethis::create_package('C:/Users/eager/OneDrive - Northeastern University/Desktop/NORTHEASTERN/Spring 2025/MMLN Package Testing/MMLN')

## 3. Update/Organize the R Files and place them in the requisite folders
## README.md will go in MMLN/
## All other .R files, with functions and definitions, will go in MMLN/R/

## 4. Add Package Metadata in DESCRIPTION

## 5. Document functions with Roxygen2
# Need to go through each .R files and add roxygen2 comments directly above
# This includes @import or @importFrom tags for functions from other packages
# When done documenting, type
# devtools::document()
# which will parse through all the roxygen blocks and update NAMESPACE and .Rd help files

## 6. Update the README.md and place it in the top-level (see step 3)
# Add the following two lines to install from GitHub (once MMLN is hosted):
# # install.packages("devtools")
# devtools::install_github("eaegerber/MMLN")
# Also add an example and a quick usage snippet

## 7. Check and Build the package
# need to make sure everything meets CRAN standards with
# devtools::check()
# If you see 0 errors, 0 warnings, 0 notes, all good. Otherwise address any issues
# You can then install the package locally with:
# devtools::install()
# library(MMLN)

## 8. Version Control and GitHub
# Create a new GitHub repo called MMLN (see step 6)
# In a terminal, initialize the repo, commit the files, and push
# git init
# git add .
# git commit -m "Initial commit"
# git branch -M main
# git remote add origin https://github.com/eaegerber/MMLN.git
# git push -u origin main
# Then, going forward, each change made to the package can be committed and pushed

## 9. (Optional) Add Tests
# can set up automated tests with a tests/testthat/ structure
# usethis::use_testthat()
# this creates a skeleton in tests/, then can create test files (e.g. test_MMLN.R)
# usethis::use_test("MMLN")
# Fill it with code to verify your functions, something like:
# test_that("MMLN runs on example data", {
#  fit <- MMLN(formula = ..., data = ...)
#  expect_true(is.list(fit))
#})
# And then rul all tests with:
# devtools::test()

## 10. Publishing and Sharing
# Once everything is on GitHub, you're essentially public
# But to submit to Cran:
# Ensure devtools::check() returns 0 errors, 0 warnings, 0 notes
# Build a submission tarball:
# devtools::build()
# This should generate: MMLN_0.0.0.9000.tar.gz
# And submit via CRAN submission portal:
