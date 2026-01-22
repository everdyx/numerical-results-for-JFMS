package main;

import java.io.*;
import java.util.ArrayList;
import java.util.List;


public class CSVProcess {
    public int getNumLines(String filepath) throws IOException {
        File file = new File(filepath);

        BufferedReader br = new BufferedReader(new FileReader(file));
        int rowCount = 0;
        while (br.readLine() != null) { // 这将读取文件的每一行，直到没有更多的行为止
            rowCount++;
        }
        br.close();

        return rowCount;
    }

    public void writeHead(ArrayList<String> header, String filepath) throws IOException {
        FileWriter writer = new FileWriter(filepath);
        // 写入表头
        for (String head:header ) {
            writer.append(head);
            writer.append(",");
        }
        writer.append("\n");
        writer.close();
    }

    public void writeALine(ArrayList<String> line, String filepath) throws IOException {
        FileWriter writer = new FileWriter(filepath,true);
        // 写入表头
        for (String a : line ) {
            writer.append(a);
            writer.append(",");
        }
        writer.append("\n");
        writer.close();
    }

    public void replaceALine(ArrayList<String> line, String filepath,int rowToRep) throws IOException {
        // 1. 读取CSV文件到一个列表（List）中
        BufferedReader br = new BufferedReader(new FileReader(filepath));
        List<String[]> rows = new ArrayList<>();
        String l;
        while ((l= br.readLine()) != null) {
            String[] row = l.split(",");
            rows.add(row);
        }
        br.close();

        // 2. 从列表中删除需要删除的行
        rows.remove(rowToRep-1);
        String[] x=new String[line.size()];
        for (int i = 0; i < x.length; i++) {
            String str=line.get(i);
            x[i]=str;
        }
        rows.add(rowToRep-1,x);

        // 3. 将更新后的列表写回到CSV文件中
        FileWriter fw = new FileWriter(filepath);
        for (String[] row : rows) {
            fw.append(String.join(",", row));
            fw.append("\n");
        }
        fw.close();
    }
}
