package algorithm

import (
	"container/heap"
	"math"
	"math/rand"
	"time"
)

//MAX_INT
const max = int64(^uint64(0) >> 1)

/*
	@note:This function serves as a mapgenerator which generating test environment
	for search algorithm.

	@input:Parameter `m` and `n` denote the row number and col number of the grid.
	Parameter `dense` denotes the blocked cell of grid. `costL` and `costH` range
	the open low boundary and close upper boundary of the process of cost generation.

	@output:`retFeasibel` returns the map of a task, and `retCost` returns the cost for an
	angent moving to the current cell from its cell neigbors.
*/
func MapGenerator(m int, n int, dense float64, costL int, costH int) (retFeasible [][]int, retCost [][]int) {
	totalCount := m * n
	//Feasible Matrix
	blockCount := int(float64(m*n) * dense)
	bucketA := make([]int, totalCount)
	for i := 0; i < blockCount; i++ {
		bucketA[i] = 1
	}
	for i := blockCount; i < totalCount; i++ {
		bucketA[i] = 0
	}

	knuthDurstenfeldShuffle(bucketA, 10)
	for i := 0; i < m; i++ {
		arr := make([]int, n)
		for j, _ := range arr {
			arr[j] = bucketA[0]
			bucketA = bucketA[1:]
		}
		retFeasible = append(retFeasible, arr)
	}
	//Cost Matrix
	bucketA = make([]int, totalCount)
	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	for i, _ := range bucketA {
		bucketA[i] = costL + r.Intn(costH-costL)
	}
	for i := 0; i < m; i++ {
		arr := make([]int, n)
		for j, _ := range arr {
			arr[j] = bucketA[0]
			bucketA = bucketA[1:]
		}
		retCost = append(retCost, arr)
	}
	return
}

/*
	@note:Here feasibleG is a two dimension slice with elements
	zeros and ones denote unblock cells and block cells respectively,
	and costG in the same type denote the current G cost in A* Algorithem.
	Furthermore, when element of feasibleG become 2, it means that the point
	have been relaxed by djkstra algorithem. It returns the least cost and
	least cost for points relaxed(dynamic programming matrix).

	@input:

	@output:
*/
func DijkstraForGrid(feasibleMap [][]int, costG [][]int, start [2]int, target [2]int) (retCost int64, step int64, tract [][2]int) {
	feasibleG := make([][]int, 0)
	for _, v := range feasibleMap {
		arr := make([]int, len(feasibleMap[0]))
		for k, _ := range v {
			arr[k] = v[k]
		}
		feasibleG = append(feasibleG, arr)
	}
	setA := make(PriorityQueue, 0)
	m := len(feasibleG)
	if m == 0 {
		return
	}
	n := len(feasibleG[0])
	dpG := make([][]int64, 0)
	//init the dp value
	for i := 0; i < m; i++ {
		arr := make([]int64, n)
		for j, _ := range arr {
			arr[j] = int64(max)
		}
		dpG = append(dpG, arr)
	}
	dpG[start[0]][start[1]] = 0
	feasibleG[target[0]][target[1]] = 0 //somewhat gurantee the feasible
	initItem := Node{
		Cost:   0,
		Id:     start[0]*n + start[1],
		Parent: nil,
	}
	heap.Push(&setA, &initItem)
	//BFS for djkstra, 8 connected
	// rowNeigbor, colNeigbor := []int{0, 1, -1}, []int{0, 1, -1}
	rowNeigbor, colNeigbor := []int{1, -1}, []int{1, -1}
	for len(setA) > 0 {
		relaxPoint := heap.Pop(&setA).(*Node)
		step++
		id := relaxPoint.Id
		nowRow, nowCol := id/n, id%n
		if feasibleG[nowRow][nowCol] == 2 {
			continue
		}
		feasibleG[nowRow][nowCol] = 2

		if nowRow == target[0] && nowCol == target[1] {
			retCost = dpG[nowRow][nowCol]
			tractCompress := extract(relaxPoint)
			tract = decodeTract(tractCompress, n)
			return
		}
		for _, v1 := range rowNeigbor {

			tempX, tempY := nowRow+v1, nowCol
			//inside the boundary & unblocked & unrelaxed
			if tempX < m && tempX >= 0 && tempY < n && tempY >= 0 && feasibleG[tempX][tempY] == 0 {
				// diag := math.Abs(float64(v1*v2))*0.4 + 1.0
				cost := dpG[nowRow][nowCol] + int64(costG[tempX][tempY])
				dpG[tempX][tempY] = cost
				item := Node{
					Cost:   cost,
					Id:     tempX*n + tempY,
					Parent: relaxPoint,
				}
				heap.Push(&setA, &item)
			}
		}

		for _, v1 := range colNeigbor {
			tempX, tempY := nowRow, nowCol+v1
			//inside the boundary & unblocked & unrelaxed
			if tempX < m && tempX >= 0 && tempY < n && tempY >= 0 && feasibleG[tempX][tempY] == 0 {
				cost := min(dpG[nowRow][nowCol]+int64(costG[tempX][tempY]), dpG[tempX][tempY])
				dpG[tempX][tempY] = cost
				item := Node{
					Cost:   cost,
					Id:     tempX*n + tempY,
					Parent: relaxPoint,
				}
				heap.Push(&setA, &item)
			}
		}
	}
	return
}

/*
	@note:A* has extra interface as a function hvalue to calculate the h-value.
	As for hvalue, it requires four parameters denote current row index, current col index,
	target row index and target col index respectively which return the h value. Other
	parameters share the same fucntion as DjkstraForGrid.

	@input:

	@output:
*/
func AstarSearch(feasibleMap [][]int, costG [][]int, start [2]int, target [2]int, hvalue func(int, int, int, int) int64) (fcost int64, step int64, tract [][2]int) {
	feasibleG := make([][]int, 0)
	for _, v := range feasibleMap {
		arr := make([]int, len(feasibleMap[0]))
		for k, _ := range v {
			arr[k] = v[k]
		}
		feasibleG = append(feasibleG, arr)
	}
	setA := make(PriorityQueue, 0)
	m := len(feasibleG)
	if m == 0 {
		return 0, 0, nil
	}
	n := len(feasibleG[0])
	dpG := make([][]int64, 0)
	//init the dp value
	for i := 0; i < m; i++ {
		arr := make([]int64, n)
		for j, _ := range arr {
			arr[j] = max
		}
		dpG = append(dpG, arr)
	}
	dpG[start[0]][start[1]] = 0
	feasibleG[target[0]][target[1]] = 0 //somewhat gurantee the feasible
	initItem := Node{
		Cost:   0,
		Id:     start[0]*n + start[1],
		Parent: nil,
		Gvalue: 0,
		Hvalue: hvalue(start[0], start[1], target[0], target[1]),
	}
	heap.Push(&setA, &initItem)
	//BFS for djkstra
	// rowNeigbor, colNeigbor := []int{0, 1, -1}, []int{0, 1, -1} //8-connected,remains to be construct
	rowNeigbor, colNeigbor := []int{1, -1}, []int{1, -1} //4-connected
	for len(setA) > 0 {
		relaxPoint := heap.Pop(&setA).(*Node)
		step++
		id := relaxPoint.Id
		nowRow, nowCol := id/n, id%n

		if nowRow == target[0] && nowCol == target[1] {
			fcost = dpG[nowRow][nowCol]
			tractCompress := extract(relaxPoint)
			tract = decodeTract(tractCompress, n)
			return
		}

		for _, v1 := range rowNeigbor {
			tempX, tempY := nowRow+v1, nowCol
			//inside the boundary & unblocked & unrelaxed
			if tempX < m && tempX >= 0 && tempY < n && tempY >= 0 && feasibleG[tempX][tempY] == 0 {
				feasibleG[tempX][tempY] = 2
				//gval := min(dpG[nowRow][nowCol]+costG[tempX][tempY], dpG[tempX][tempY])
				gval := dpG[nowRow][nowCol] + int64(float64(costG[tempX][tempY]))
				hval := hvalue(tempX, tempY, target[0], target[1])
				dpG[tempX][tempY] = gval
				item := Node{
					Cost:   gval + hval,
					Id:     tempX*n + tempY,
					Parent: relaxPoint,
					Gvalue: gval,
					Hvalue: hval,
				}
				heap.Push(&setA, &item)

			}
		}
		for _, v1 := range colNeigbor {
			tempX, tempY := nowRow, nowCol+v1
			//inside the boundary & unblocked & unrelaxed
			if tempX < m && tempX >= 0 && tempY < n && tempY >= 0 && feasibleG[tempX][tempY] == 0 {
				feasibleG[tempX][tempY] = 2
				//gval := min(dpG[nowRow][nowCol]+costG[tempX][tempY], dpG[tempX][tempY])
				gval := dpG[nowRow][nowCol] + int64(costG[tempX][tempY])
				hval := hvalue(tempX, tempY, target[0], target[1])
				dpG[tempX][tempY] = gval
				item := Node{
					Cost:   gval + hval,
					Id:     tempX*n + tempY,
					Parent: relaxPoint,
					Gvalue: gval,
					Hvalue: hval,
				}
				heap.Push(&setA, &item)
			}
		}

	}
	return
}

/*
	@note: This function use concurency method to search the route in forward and backward direction,
	it's useful in task such as agents to find each other.

	@input:

	@output:
*/

//Alternative in function A*, the F value evaluation method
func HalmintanDistance(currentX, currentY, targetX, targetY int) int64 {
	abs1 := currentX - targetX
	if abs1 < 0 {
		abs1 = -abs1
	}
	abs2 := currentY - targetY
	if abs2 < 0 {
		abs2 = -abs2
	}
	return int64(abs1 + abs2)
}

//Alternative in function A*, the F value evaluation method
func EulerDistance(currentX, currentY, targetX, targetY int) int64 {
	abs1 := currentX - targetY
	abs2 := targetX - currentY
	return int64(math.Sqrt(float64(abs1*abs1 + abs2*abs2)))
}

//Alternative in function A*, the F value evaluation method
func ChebyshevDistance(currentX, currentY, targetX, targetY int) int64 {
	abs1 := currentX - targetX
	if abs1 < 0 {
		abs1 = -abs1
	}
	abs2 := currentY - targetY
	if abs2 < 0 {
		abs2 = -abs2
	}
	return int64(float64(abs1+abs2) + 1.414213*float64(min(int64(abs1), int64(abs2))))
}

//Shuffling data for random tasks
func knuthDurstenfeldShuffle(list []int, additionNum int) {
	l := len(list)
	if l <= 1 {
		return
	}
	r := rand.New(rand.NewSource(time.Now().UnixNano() + int64(additionNum<<1)))
	for l > 0 {
		randIdx := r.Intn(l)
		list[l-1], list[randIdx] = list[randIdx], list[l-1]
		l--
	}
}

//util
func min(i, j int64) int64 {
	if i < j {
		return i
	} else {
		return j
	}
}

//imply to extract the route after calculate
func extract(target *Node) []int {
	ret := make([]int, 0)
	for target != nil {
		ret = append(ret, target.Id)
		target = target.Parent
	}
	return ret
}

//imply to decode the compress route of algorithm
func decodeTract(tract []int, n int) (ret [][2]int) {
	for _, v := range tract {
		temp := [2]int{v / n, v % n}
		ret = append(ret, temp)
	}
	return
}
