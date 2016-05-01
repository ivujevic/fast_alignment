## Building

### Build requirements

One needs `GNU make`, `Cmake 3+`, `Boost` libraries. If you don't know if you have them or
not just type `make packages` and it will prompt you for install.

Built it with Boost 1.60 but it should work with older versions if interface
didn't change.

### How to build?

Build by writing `make` in your command line. If everything went without errors
executables should be present in the directory called `build`. There is executable file
tachyonn

## Running

There are three executables. If you want to use this app on client-server arhitecture use `tachyon_server`.
It dump base in memory, waiting for queries and store outputs to given path. For this, you also need other executable.
In order to run this procedure you need to start `run_query` with appropriate parameters. If you want to use this app
as stand-alone application you will need third executable, `tachyon`.

