package main

import algorithm "basic/models/algorithm/searchmethod"

func main() {
	//tract调试
	// f, _ := os.Create("trace.out")
	// defer f.Close()
	// trace.Start(f)
	algorithm.Compare(50, 50, 0.2, 10, 11, [2]int{1, 1}, [2]int{49, 49}, 0)
	// trace.Stop()
}
