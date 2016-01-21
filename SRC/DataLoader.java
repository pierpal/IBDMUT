/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
 */

package IBDMUT;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.BitSet;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Pier Palamara
 */
public class DataLoader implements Callable<Dataset> {

    private final String pedFile;
    private final boolean saveBin;
    private final boolean loadBin;
    private final double posteriorFrom;
    private final double posteriorTo;
    private final boolean storePosterior;
    private final boolean usePosteriors;
    private final String saveBinSuffix;
    private final String loadBinSuffix;

    public DataLoader(String loadBinSuffix, String saveBinSuffix, String pedFile, boolean saveBin, boolean loadBin, boolean usePosteriors, double posteriorFrom, double posteriorTo, boolean storePosterior) {
        this.loadBin = loadBin;
        this.pedFile = pedFile;
        this.saveBin = saveBin;
        this.posteriorFrom = posteriorFrom;
        this.posteriorTo = posteriorTo;
        this.storePosterior = storePosterior;
        this.usePosteriors = usePosteriors;
        this.loadBinSuffix = loadBinSuffix;
        this.saveBinSuffix = saveBinSuffix;
    }

    public Dataset call() throws FileNotFoundException, IOException {
        Dataset data;
        if (!loadBin) {
            data = new Dataset();
            try {
                data.readPedFile(pedFile);
            } catch (Exception ex) {
                Logger.getLogger(PedMatchProcessor.class.getName()).log(Level.SEVERE, null, ex);
                System.exit(1);
            }
            File f = new File(pedFile + ".frq");
            if (!f.exists()) {
                Tools.exit("Frequency file not " + f + "found");
            } else {
                data.addFreqFile(pedFile + ".frq");
                data.setHaveFreq(true);
            }

            if (usePosteriors) {
                String posteriorsFile = pedFile + ".post.gz";
                try {
                    data.readPosteriorsFile(posteriorsFile, posteriorFrom, posteriorTo, storePosterior);
                } catch (Exception ex) {
                    Logger.getLogger(PedMatchProcessor.class.getName()).log(Level.SEVERE, null, ex);
                    System.exit(1);
                }
            }
            if (saveBin) {
                String binFile = pedFile + saveBinSuffix + ".bin.gz";
                Tools.printVerboseProgressLevel2("Saving data to " + binFile);
                Kryo kryo = new Kryo();
                kryo.register(BitSet.class, new BitSetSerializer());
                OutputStream outputStream = new GZIPOutputStream(new FileOutputStream(binFile));
                Output output = new Output(outputStream);
                kryo.writeObject(output, data);
                output.close();
            }
        } else {
            String binFile = pedFile + loadBinSuffix + ".bin.gz";
            Tools.printProgress("Loading data from " + binFile);
            Kryo kryo = new Kryo();
            kryo.register(BitSet.class, new BitSetSerializer());
            GZIPInputStream inputStream = new GZIPInputStream(new FileInputStream(binFile));
            Input input = new Input(inputStream);
            data = kryo.readObject(input, Dataset.class);
            input.close();
        }
        return data;
    }

}
