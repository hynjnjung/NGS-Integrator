package nih.nhlbi.esbl.ngs;

import java.util.Comparator;
import java.util.PriorityQueue;

class MedianFinder {

    PriorityQueue<Double> minHeap = null;
    PriorityQueue<Double> maxHeap = null;

    public MedianFinder() {
        minHeap = new PriorityQueue<>();
        maxHeap = new PriorityQueue<>(Comparator.reverseOrder());
    }

    public void addNum(Double num) {
        minHeap.offer(num);
        maxHeap.offer(minHeap.poll());

        if (minHeap.size() < maxHeap.size()) {
            minHeap.offer(maxHeap.poll());
        }
    }

    public void removeNum(Double num) {
        if (minHeap.contains(num)) {
            minHeap.remove(num);

            //Rebalancing heaps
            if (minHeap.size() < maxHeap.size()) {
                minHeap.offer(maxHeap.poll());
            }

        } else {
            maxHeap.remove(num);

            //Rebalancing heaps
            if (minHeap.size() > maxHeap.size() + 1) {
                maxHeap.offer(minHeap.poll());
            }
        }

    }

    public double findMedian() {
        if (minHeap.size() > maxHeap.size()) {
            return minHeap.peek();
        } else {
            return (minHeap.peek() + maxHeap.peek()) / 2.0;
        }
    }
}
