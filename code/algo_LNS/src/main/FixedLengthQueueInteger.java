package main;

import java.util.ArrayDeque;

public class FixedLengthQueueInteger {
    private final int MAX_SIZE = 2; // 固定长度
    private ArrayDeque<Integer> queue = new ArrayDeque<>(MAX_SIZE);

    public void enqueue(int value) {
        if (queue.size() == MAX_SIZE) {
            // 如果队列已满，移除头部元素
            queue.removeFirst();
        }
        // 添加新元素到队尾
        queue.addLast(value);
    }

    public int dequeue() {
        if (queue.isEmpty()) {
            throw new IllegalStateException("Queue is empty");
        }
        // 移除并返回队首元素
        return queue.removeFirst();
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