build:
    if test -d deploy; then rm -rf deploy; fi
    mkdir deploy
    wasm-pack build --target web --out-dir ./static
    cp ./static/*.js deploy
    cp ./static/*.wasm deploy
    cp index.html deploy
zip:
    zip -r -j deploy.zip deploy
