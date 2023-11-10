CC= gcc
CFLAGS= -g -O2 -Wall
POSTFLAGS= -lm 

kekule: main.o input.o concat.o fold.o signature.o avl.o index.o chrono.o output.o 
	$(CC) $(CFLAGS) -o kekule main.o input.o concat.o fold.o signature.o avl.o index.o chrono.o output.o $(POSTFLAGS)

main.o: main.c structure.h input.h concat.h avl.h
	$(CC) $(CFLAGS) -c main.c 

input.o: input.c input.h structure.h
	$(CC) $(CFLAGS) -c input.c

concat.o: concat.c concat.h structure.h input.h fold.h avl.h
	$(CC) $(CFLAGS) -c concat.c 

fold.o: fold.c fold.h structure.h concat.h signature.h avl.h
	$(CC) $(CFLAGS) -c fold.c

signature.o: signature.c signature.h structure.h
	$(CC) $(CFLAGS) -c signature.c 

avl.o: avl.c avl.h signature.h structure.h signature.h
	$(CC) $(CFLAGS) -c avl.c 

index.o: index.c index.h
	$(CC) $(CFLAGS) -c index.c 

chrono.o: chrono.c chrono.h
	$(CC) $(CFLAGS) -c chrono.c 

output.o: output.c output.h
	$(CC) $(CFLAGS) -c output.c

clean:
	rm -f *.o
	rm -f kekule
	rm -rf Kekule kekule.tgz

mrproper: clean
	rm -f callgrind.out.*

debugmeta: clean kekule
	echo "DEBUGMETA"
	valgrind ./kekule Examples_2/testYI.mot 10 Examples_2/testYI.conc

meta: clean kekule
	echo "META"
	./kekule Examples_2/testYI.mot 20 Examples_2/testYIsimple.conc
	./kekule Examples_2/testXI.mot 18 Examples_2/testXI.conc

YI: clean kekule
	./kekule Examples_2/testYI.mot 10
	./kekule Examples_2/testYI.mot 15
	./kekule Examples_2/testYI.mot 20
	./kekule Examples_2/testYI.mot 25	

YVV: clean kekule
	./kekule Examples_2/testYVV.mot 8
	./kekule Examples_2/testYVV.mot 16
	./kekule Examples_2/testYVV.mot 24

XVI:
	./kekule Examples_2/testXVI.mot 4
	./kekule Examples_2/testXVI.mot 8
	./kekule Examples_2/testXVI.mot 12
	./kekule Examples_2/testXVI.mot 16
	./kekule Examples_2/testXVI.mot 20

IVV: clean kekule
	./kekule Examples_2/testIVV.mot 6
	./kekule Examples_2/testIVV.mot 9
	./kekule Examples_2/testIVV.mot 12
	./kekule Examples_2/testIVV.mot 15
	./kekule Examples_2/testIVV.mot 18
	./kekule Examples_2/testIVV.mot 24

JYV: clean kekule
	./kekule Examples_2/testJYV.mot 4
	./kekule Examples_2/testJYV.mot 8
	./kekule Examples_2/testJYV.mot 12
	./kekule Examples_2/testJYV.mot 16
	./kekule Examples_2/testJYV.mot 20
	
YVI: clean kekule
	./kekule Examples_2/testYVI.mot 13

Y4: clean kekule
	./kekule Examples_2/testY4.mot 4
	./kekule Examples_2/testY4.mot 6
	./kekule Examples_2/testY4.mot 8
	./kekule Examples_2/testY4.mot 10

XY: 
	./kekule Examples_2/testXY.mot 7
	./kekule Examples_2/testXY.mot 14
	./kekule Examples_2/testXY.mot 21

IVK:
	./kekule Examples_2/testIVK.mot 6
	./kekule Examples_2/testIVK.mot 12
	./kekule Examples_2/testIVK.mot 18


XI: clean kekule
	./kekule Examples_2/testXI.mot 6
	./kekule Examples_2/testXI.mot 9
	./kekule Examples_2/testXI.mot 12
	./kekule Examples_2/testXI.mot 15
	./kekule Examples_2/testXI.mot 18

JVV: clean kekule
	./kekule Examples_2/testJVV.mot 6
	./kekule Examples_2/testJVV.mot 9	
	./kekule Examples_2/testJVV.mot 12
	./kekule Examples_2/testJVV.mot 15
		

valmeta: clean kekule
	valgrind --leak-check=full --show-leak-kinds=all ./kekule Examples_2/testXI-modif.mot 3 Examples_2/testXI-modif.trans

debug: clean kekule
	valgrind --leak-check=full --show-leak-kinds=all ./kekule Examples_2/testIVV.mot 9

profile: clean kekule
	valgrind --tool=callgrind ./kekule Examples_2/testIVV.mot 12
	kcachegrind

git: clean
	git pull


tar: clean
	mkdir Kekule
	cp *.c *.h Kekule
	cp Makefile Kekule
	cp -r Examples_2 Kekule
	tar cvfz kekule.tgz Kekule
	rm -rf Kekule
	ls -sh kekule.tgz
