library(usethis)
use_data_raw()
use_package("mice")
use_package("truncdist")
use_package("nlshrink")
use_package("MASS")
use_package("glmnet")
use_package("missMDA")
use_package("SLOPE")
use_package("stats")

use_rcpp()
#Generate the necessary modifications to your NAMESPACE by documenting them with Ctrl/Cmd + Shift + D.

#Click Build & Reload in the build pane, or press Ctrl/Cmd + Shift + B.


compileAttributes(verbose=TRUE)
use_rcpp_armadillo()
library(tools)
package_native_routine_registration_skeleton(".")


