package main;

import java.util.ArrayDeque;

public class FixedLengthQueueDouble {
    private final int MAX_SIZE = 2; // 固定长度
    private ArrayDeque<Double> queue = new ArrayDeque<>(MAX_SIZE);

    public void enqueue(double value) {
        if (queue.size() == MAX_SIZE) {
            // 如果队列已满，移除头部元素
            queue.removeFirst();
        }
        // 添加新元素到队尾
        queue.addLast(value);
    }

    public double dequeue() {
        if (queue.isEmpty()) {
            throw new IllegalStateException("Queue is empty");
        }
        // 移除并返回队首元素
        return queue.removeFirst();
    }

    public double getMin() {
        // 返回队首元素
        return Math.min(queue.getFirst(),queue.getLast());
    }
    public double getLast() {
        // 返回队wei元素
        return queue.getLast();
    }

    public boolean isEmpty() {
        return queue.isEmpty();
    }

    // public static void main(String[] args) {
    //     FixedLengthQueue fifoQueue = new FixedLengthQueue();
    //     for (int i = 1; i <= 15; i++) {
    //         fifoQueue.enqueue(i);
    //         System.out.println("After adding " + i + ": " + fifoQueue.queue);
    //     }
    // }
}