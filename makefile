# コンパイラとフラグの設定
FC = gfortran
FFLAGS = -O2 -Wall

# ターゲット名とオブジェクトファイル
TARGET = solve_electric_potential
OBJS = main.o variables_module.o solve_electric_potential.o

# デフォルトターゲット
all: $(TARGET)

# ターゲットのリンク
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS)

# モジュールのコンパイル
variables_module.o: variables_module.f90
	$(FC) $(FFLAGS) -c variables_module.f90

# メインプログラムのコンパイル
solve_electric_potential.o: solve_electric_potential.f90 variables_module.o
	$(FC) $(FFLAGS) -c solve_electric_potential.f90

# メインプログラムのコンパイル
main.o: main.f90 variables_module.o
	$(FC) $(FFLAGS) -c main.f90

# クリーンアップターゲット
clean:
	rm -f $(OBJS) $(TARGET)

# 実行ターゲット
run: $(TARGET)
	./$(TARGET)
