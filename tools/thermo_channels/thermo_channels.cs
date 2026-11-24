// List the non-MS detectors stored in a Thermo .raw file (analog, UV, PDA, MS analog) and export
// each detector's traces as CSV (time in minutes, one column per channel).
// Uses Thermo's RawFileReader library, which ships with ThermoRawFileParser.
//
// Build:  mcs -r:<dir>/ThermoFisher.CommonCore.Data.dll -r:<dir>/ThermoFisher.CommonCore.RawFileReader.dll
//             thermo_channels.cs -out:thermo_channels.exe
// Run:    MONO_PATH=<dir> mono thermo_channels.exe run.raw [output_folder]
using System;
using System.IO;
using System.Linq;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

class ThermoChannels
{
    static int Main(string[] args)
    {
        if (args.Length < 1) { Console.Error.WriteLine("usage: thermo_channels.exe run.raw [out]"); return 2; }
        string outDir = args.Length > 1 ? args[1] : null;
        if (outDir != null) Directory.CreateDirectory(outDir);
        IRawDataPlus raw = RawFileReaderAdapter.FileFactory(args[0]);
        if (raw.IsError || !raw.IsOpen) { Console.Error.WriteLine("cannot open " + args[0]); return 1; }
        raw.IncludeReferenceAndExceptionData = false;
        Console.WriteLine("file: " + Path.GetFileName(args[0]));
        Device[] kinds = { Device.MS, Device.Analog, Device.UV, Device.Pda, Device.MSAnalog, Device.Other };
        foreach (Device kind in kinds)
        {
            int n = raw.GetInstrumentCountOfType(kind);
            Console.WriteLine($"{kind}: {n} device(s)");
            if (kind == Device.MS) continue;
            for (int k = 1; k <= n; k++)
            {
                raw.SelectInstrument(kind, k);
                var info = raw.GetInstrumentData();
                var hdr = raw.RunHeaderEx;
                string[] labels = (info.ChannelLabels ?? new string[0]).Where(s => !string.IsNullOrEmpty(s)).ToArray();
                Console.WriteLine($"  {kind} {k}: name '{info.Name}', model '{info.Model}', units {info.Units}, " +
                                  $"channels [{string.Join(", ", labels)}], scans {hdr.FirstSpectrum}-{hdr.LastSpectrum}, " +
                                  $"time {hdr.StartTime:F2}-{hdr.EndTime:F2} min");
                if (outDir == null || hdr.LastSpectrum < hdr.FirstSpectrum) continue;
                string path = Path.Combine(outDir, $"{Path.GetFileNameWithoutExtension(args[0])}_{kind}{k}.csv");
                using (var w = new StreamWriter(path))
                {
                    int width = -1;
                    for (int s = hdr.FirstSpectrum; s <= hdr.LastSpectrum; s++)
                    {
                        double t = raw.RetentionTimeFromScanNumber(s);
                        var seg = raw.GetSegmentedScanFromScanNumber(s, null);
                        double[] v = seg.Intensities ?? new double[0];
                        if (width < 0)
                        {
                            width = v.Length;
                            var head = Enumerable.Range(0, width).Select(i => i < labels.Length ? labels[i] : $"channel{i + 1}");
                            w.WriteLine("time_min," + string.Join(",", head));
                        }
                        w.WriteLine(t.ToString("R") + "," + string.Join(",", v.Take(width).Select(x => x.ToString("R"))));
                    }
                }
                Console.WriteLine($"    written {path}");
            }
        }
        raw.Dispose();
        return 0;
    }
}
