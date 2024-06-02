# This sample should be run after make libHEAAN.a file in lib folder

# All Target
all: clean HEAAN MY_TEST run

# Tool invocations
HEAAN:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	# g++ -L/usr/local/lib -o "HEAAN" ../lib/libHEAAN.a run.cpp -I/../src -lgmp -lntl -lm -O2 -std=c++11  这行代码是作者库里写的，有问题，用下一行master里的代码就可以生成
	g++ run.cpp -o HEAAN -L/usr/local/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	@echo 'Finished building target: $@'
	@echo ' '

MY_TEST:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ my_test.cpp -o MY_TEST -L/usr/local/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	@echo 'Finished building target: $@'
	@echo ' '	

# Other Targets
clean:
	rm -rf HEAAN MY_TEST

run:
	# ./HEAAN
	./MY_TEST