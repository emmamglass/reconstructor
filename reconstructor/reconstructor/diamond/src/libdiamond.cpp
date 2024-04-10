#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdexcept>
#include "./run/main.h"
#include "./basic/const.h"
#include "./util/io/exceptions.h"

static PyObject* method_run_blastp(PyObject* self, PyObject* args)
{
    try {
        int size = PyTuple_Size(args);
        int argc = size + 2;
        const char** argv = new const char*[argc];
        argv[0] = "diamond";
        argv[1] = "blastp";
        for (int i = 0; i < size; i++) {
            PyObject* item = PyTuple_GetItem(args, i);
            if (!PyArg_Parse(item, "s", &(argv[i + 2]))) {
                delete[] argv;
                return NULL;
            }
        }
        int status = main(argc, argv);
        delete[] argv;
        return Py_BuildValue("i", status);
	}
	catch(const std::bad_alloc &e) {
        PyErr_SetString(PyExc_MemoryError, e.what());
        return NULL;
    }
    catch (const FileOpenException& e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
	}
    catch (const File_read_exception& e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }
    catch (const File_write_exception& e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }
    catch (const EndOfStream& e) {
        PyErr_SetString(PyExc_EOFError, e.what());
        return NULL;
    }
    catch (const StreamReadException& e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }
    catch(const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
    catch(...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown error occurred in extension");
        return NULL;
    }
};

static PyObject* method_version(PyObject *self, PyObject*args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }
    char* version = (char *)Const::version_string;
    return Py_BuildValue("s", version);
}

static PyMethodDef libdiamond_methods[] = {
    {
        "run_blastp", 
        method_run_blastp, 
        METH_VARARGS, 
        "Run diamond by its command options. For option details, just pass a 'help' argument"
    },
    {
        "version", 
        method_version, 
        METH_VARARGS, 
        "Return the version of diamond."
    },
    {
        NULL, 
        NULL, 
        0, 
        NULL}
};

static struct PyModuleDef libdiamond_module {
    PyModuleDef_HEAD_INIT,
    "libdiamond",
    "Python wrapper for Diamond supporting Reconstructor",
    -1,
    libdiamond_methods
};

PyMODINIT_FUNC PyInit_libdiamond(void)
{
    return PyModule_Create(&libdiamond_module);
}