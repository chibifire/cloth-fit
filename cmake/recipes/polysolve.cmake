# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT

if(TARGET polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve'")

include(CPM)
CPMAddPackage("gh:Huangzizhou/polysolve#ee553dc6a330ef70bcc6ba2f05653fecd201c05f")
