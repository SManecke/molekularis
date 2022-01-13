TEXMFOUTPUT:=build
TEXSOURCES=$(wildcard *.tex)
all: $(patsubst %.tex,build/%.pdf,$(TEXSOURCES))
	mkdir -p build
	clang -O3 -Wall -Wextra main.c picosat.c -lm -o build/molecularis
	clang -g -O0 -Wall -Werror gui.c  -lraylib -ldl -lX11 -lglfw -lpthread -lm -o build/gui

	clang -g3 -O0 -Wall -Werror wasm.c -target wasm32 -nostdlib           \
		-fvisibility=hidden -fno-builtin -fno-exceptions              \
		-fno-threadsafe-statics -Wl,--no-entry                        \
		-Wl,--allow-undefined-file=wasm.syms,,--initial-memory=131072 \
		-Wl,--export=init,--export=update,--export=set_value          \
		-Wl,--export=get_value,--export=undo,--export=redo -o build/out.wasm
	cp build/out.wasm html/molekularis.wasm

build/%.pdf: %.tex
	mkdir -p build
	latexmk -bibtex -pdf -jobname=build/$(patsubst %.tex,%,$<)            \
		-pdflatex="pdflatex -interaction=nonstopmode -halt-on-error -file-line-error" -use-make $<
