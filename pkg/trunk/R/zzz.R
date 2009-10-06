.First.lib <- function(lib, pkg) {
	library.dynam("dairy", pkg, lib)
	.Call("DAIRY_init")
}
