# EDF-Reader

```bash
  mkdir build
  cd build
  cmake ../source
  cmake --build .
```


### Compilación para WebAssembly

Para compilar el proyecto para WebAssembly, sigue estos pasos:

1. Navega al directorio `source`:

    ```bash
    cd source
    ```

2. Ejecuta el siguiente comando de Emscripten para compilar el código fuente:

    ```bash
    emcc edfwasm.cpp edf_bindings.cpp signaldata.cpp -o edf.js -s WASM=1 -s FORCE_FILESYSTEM=1 -s ALLOW_MEMORY_GROWTH=1 -s "EXPORTED_RUNTIME_METHODS=['ccall', 'cwrap']" --bind -std=c++20 -I ../Eigen
    ```