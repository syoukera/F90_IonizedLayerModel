# コンパイラとフラグの設定
FC = gfortran
FFLAGS = -O2 -Wall -Jsrc
# FFLAGS = -Wall -Jsrc -fcheck=bounds -g -O0

# ターゲット名とオブジェクトファイル
TARGET = ionized_layer
OBJS = src/main.o src/variables_module.o src/solve_poisson_equation.o src/solve_ion_conservation.o

# デフォルトターゲット
all: $(TARGET)

# ターゲットのリンク
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS)

# モジュールのコンパイル
src/variables_module.o: src/variables_module.f90
	$(FC) $(FFLAGS) -c src/variables_module.f90 -o src/variables_module.o

# メインプログラムのコンパイル
src/solve_poisson_equation.o: src/solve_poisson_equation.f90 src/variables_module.o
	$(FC) $(FFLAGS) -c src/solve_poisson_equation.f90 -o src/solve_poisson_equation.o

# メインプログラムのコンパイル
src/main.o: src/main.f90 src/variables_module.o
	$(FC) $(FFLAGS) -c src/main.f90 -o src/main.o

# メインプログラムのコンパイル
src/solve_ion_conservation.o: src/solve_ion_conservation.f90 src/variables_module.o
	$(FC) $(FFLAGS) -c src/solve_ion_conservation.f90 -o src/solve_ion_conservation.o

# クリーンアップターゲット
clean:
	rm -f $(OBJS) $(TARGET) src/*.mod