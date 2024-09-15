#!/usr/bin/env python
#  Copyright 2023
#  Center for Global Discovery and Conservation Science, Arizona State University
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#
# Original authors:
#   Kathy Heidebrecht
#   Greg Asner, gasner AT asu.edu
#
# Translated to Python by:
#   Nicholas Vaughn, nickvaughn AT asu.edu
#   Elahe Jamalinia, ejamalin AT asu.edu
#
# AutoMCU
# Code is provided to Planet, PBC as part of the CarbonMapper Land and Ocean
# Program. This methodology was developed in this current form by former and
# current members of the Asner Lab at GDCS. Please give proper attribution
# when using this code for publication:
#
# Asner, G. P., and K. B. Heidebrecht. 2002. Spectral unmixing of vegetation,
# soil and dry carbon cover in arid regions: Comparing multispectral and
# hyperspectral observations. IJRS 23:3939–3958.

import os
import argparse
from automcu import AutoMCU


def main():
    """
    Required and optional arguments to be included in the CLI, while running the code.
    The arguments are defined in the main function. The arguments are parsed in the __main__.py file, or used in the run function in the automcu.py file.
    """
    # 创建一个ArgumentParser对象，用于解析命令行参数
    parser = argparse.ArgumentParser(description="Run AutoMCU on an image")

    # 添加一个命令行参数，用于控制是否输出详细信息
    # 添加一个参数，用于输出详细的信息
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")


    # 添加一个参数，用于强制覆盖已存在的输出文件
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Force overwrite of existing output file",
    )


    # 添加一个参数，用于从JSON格式的配置文件中获取设置，而不是指定命令行参数
    parser.add_argument(
        "--config",
        "-c",
        default=None,
        help="Get settings from JSON-formatted config file,"
        " instead of needing to specify command line args",
    )

    # 添加一个参数，用于指定输入ENVI hdr的名称，如果它不是标准路径.hdr或splitext(path).hdr
    parser.add_argument(
        "--input_hdr",
        default=None,
        help="Specify name of input ENVI hdr in case it is not"
        " the standard path.hdr or splitext(path).hdr",
    )

    # 添加一个参数，用于指定用于解混的波段范围，该范围由两个浮点数组成，例如(600.0,750)，这些波段的
    # 数据将被合并成一个数组进行拟合。必须指定至少一个wl窗口。wl_range和band_range不能同时指定
    parser.add_argument(
        "--wl_range",
        "-w",
        default=[],
        action="append",
        help="A 2-element list or tuple of floats, i.e."
        " (600.0,750), that will be used to identify"
        " the band ranges used for unmixing."
        " Data from the bands specified by"
        " these wavelengths will be lumped into a single"
        " array for fitting. A minimum of one wl window"
        " must be specified. Both wl_range and band_range"
        " cannot be specified",
    )

    # 添加一个参数，用于指定用于解混的波段范围，该范围由两个整数组成，例如(65,72)，这些波段的
    # 数据将被合并成一个数组进行拟合。必须指定至少一个波段窗口。wl_range和band_range不能同时指定
    parser.add_argument(
        "--band_range",
        "-b",
        default=[],
        action="append",
        help="A 2-element list or tuple of integers, i.e."
        " (65,72), that will be used to identify"
        " the band ranges used for unmixing."
        " Band numbers are 0-based, so first band is 0."
        " Data from the bands specified by"
        " these windows will be lumped into a single"
        " array for fitting. A minimum of one band window"
        " must be specified. Both wl_range and band_range"
        " cannot be specified",
    )

    # 添加一个参数，用于从每个波段范围的反射光谱中减去第一个窗口值，以帮助缓解亮度问题。
    # （解混数组中的第一个值对于每个窗口的起始值变为0.0）。默认情况下这是开启的。使用此标志来禁用。
    parser.add_argument(
        "--notie",
        dest="tied",
        action="store_false",
        help="Reflectance spectra from each band range are"
        " subtracted from the first window value to"
        " help relieve brightness issues. (First value"
        " in unmixing array becomes 0.0 for each"
        " window's start). By default this is on. Use this"
        " flag to disable.",
    )

    # 添加一个参数，用于将每个波段范围的反射光谱除以第一个窗口值，以帮助缓解亮度问题。
    # （解混数组中的第一个值对于每个窗口的起始值变为1.0）。默认情况下这是关闭的。
    parser.add_argument(
        "--divide",
        action="store_true",
        help="Reflectance spectra from each band range are"
        " divided by the first window value to"
        " help relieve brightness issues. (First value"
        " in unmixing array becomes 1.0 for each"
        " window's start). By default this is off.",
    )

    # 添加一个参数，用于添加一个阴影的成分
    parser.add_argument(
        "--shade", action="store_true", help="Add an endmember for shade"
    )

    # 添加一个参数，用于在解混数组的末尾添加一个值为1.0的元素，并帮助将成分分数的总和
    # 接近1.0
    parser.add_argument(
        "--sum_to_one",
        "-1",
        action="store_true",
        help="Add an element with value 1.0 to the end of the"
        " unmixing array and the input refl to help"
        " coerce sums of the endmember fractions to be"
        " near 1.0.",
    )

    # 添加一个参数，用于指定每个像素的MC运行次数。默认为30
    parser.add_argument(
        "--iterations",
        "--count",
        type=int,
        default=30,
        help="How many MC runs to do for each pixel. " " Default 30",
    )

    # 添加一个参数，用于指定要解混的块数。默认为所有
    # 一个由两个整数组成的2元素列表或元组，例如(100,200)，用于标识用于解混的块范围

    ##define a argument having a range for a parameter --num_blocks
    parser.add_argument(
        "--num_blocks",
        default=[],
        action="append",
        help="How many blocks to unmix. "
        " Default all"
        "A 2-element list or tuple of int, i.e."
        " (100,200), that will be used to identify"
        " the blocks ranges used for unmixing.",
    )

    # 添加一个参数，用于将解混结果缩放，以使其适应输出指定的数据类型
    parser.add_argument(
        "--outscale",
        type=float,
        default=1000,
        help="Scale value applied to unmixing results to make"
        " them fit into data type specified for output.",
    )

    # 添加一个参数，用于将无数据值写入元数据。如果默认为None，输出中的缺失数据将为0，
    # 但这不会记录在元数据中
    parser.add_argument(
        "--nodata",
        type=float,
        default=None,
        help="Output no data value written into metadata. If"
        " default None, missing data in output will be 0,"
        " but this will not be recorded in metadata.",
    )

    # 添加参数，用于将输入数据乘以指定的比例，以获得0到1之间的反射率
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Apply this scale to input data to get" " reflectance scaled from 0 to 1",
    )

    # 添加参数，用于指定输出地图数据类型的Rasterio规范
    parser.add_argument(
        "--dtype",
        "--of",
        default="int16",
        help="Rasterio-usable specification of output map" " data type.",
    )

    # 添加参数，用于指定输出数据格式的GDAL短名称
    parser.add_argument(
        "--otype",
        "--ot",
        default="GTiff",
        help="GDAL shortname for output data format.",
    )

    # 添加参数，用于指定GDAL数据格式写入选项
    parser.add_argument(
        "--co",
        "-o",
        default=[],
        action="append",
        help="GDAL data format writing options.",
    )

    # 添加参数，用于假设em csvs的第二列是fwhm值。默认情况下，使用平均波段间距
    parser.add_argument(
        "--emfwhm",
        action="store_true",
        help="Assume second column of em csvs are fwhm values."
        " By default, average inter-band spacing is used",
    )

    # 添加参数，用于假设em csvs具有与输入反射率数据相同的计数和波长（将被检查）
    parser.add_argument(
        "--nointerp",
        action="store_true",
        help="Assume em csvs have same count and wavelengths"
        " as input reflectance data (will be checked)",
    )

    # 添加参数，用于重新创建原始C代码中的一个错误，该错误将输入反射率数据移动了一个索引/波段
    parser.add_argument(
        "--match_c",
        action="store_true",
        help="Recreate a bug in original C code that shifted"
        " input reflectance data by one index/band",
    )

    # 添加参数，用于为给定CSV文件定义的类提供名称。可以多次指定或使用逗号分隔的值
    parser.add_argument(
        "--names",
        default=[],
        action="append",
        metavar="STR[,STR]",
        help="Supply names for the classes defined by the"
        " given CSV files. Can be specified multiple"
        " times or using comma separated values",
    )

    # 添加参数，用于指定ENVI格式的反射率图像
    parser.add_argument("input", help="ENVI-formatted reflectance image")

    # 添加参数，用于指定AutoMCU结果地图文件名
    parser.add_argument("output", help="AutoMCU result map filename")

    # 添加参数，用于指定给定端成员类的样本列表的CSV文件。CSV应该按照样本是列，波段是行的方式排列，第一列包含波段中心波长（以纳米为单位）。
    # 至少必须指定一个类。如果文件名后面紧跟着一个:和整数n，则只读取前n个样本
    parser.add_argument(
        "emlist",
        nargs="+",
        metavar="CSV<:INT>",
        help="CSV file with list of samples for a given"
        " endmember class. CSV should should be arranged"
        " such that samples are columns, and bands are"
        " rows, with the first column containing the band"
        " center wavelengths in nanometers. At least one"
        " class must be specified. If a : and integer n"
        " immediately follow the filename, then only the"
        " first n samples will be read",
    )

    # 解析命令行参数，检查命令行参数中是否包含verbose标志，如果包含，则打印出所有命令行参数及其对应的值
    args = parser.parse_args()
    if args.verbose:
        print("Args:")
        for k, v in vars(args).items():
            print(f"{k:12s}: {v}")

    #解析一个名为 args.emlist 的列表，该列表包含了一些字符串，每个字符串可能包含一个冒号分隔的值。
    # 代码的主要目的是从这些字符串中提取两个列表：em_csvs 和 em_counts
    em_csvs = []
    em_counts = []
    for st in args.emlist:
        # 使用 strip() 方法去除字符串两端的空白字符，然后使用 split(":") 方法按照冒号分割字符串，得到一个列表 parts
        parts = st.strip().split(":")
        if len(parts) > 1:
            em_csvs.append(parts[0])
            if len(parts[1]) > 0:
                em_counts.append(int(parts[1]))
            else:
                em_counts.append(0)
        else:
            em_csvs.append(parts[0])
            em_counts.append(0)

    em_names = []
    # 首先检查 args.names 是否存在
    if args.names:
        # 如果 args.names 存在，遍历列表中的每个元素 n
        for n in args.names:
            # 对每个元素 n 使用 strip() 方法去除字符串两端的空白字符，然后使用 split(",") 方法按照逗号分割字符串，得到一个子列表 parts
            parts = n.strip().split(",")
            # 使用 extend() 方法将 parts 列表中的所有元素添加到 em_names 列表的末尾
            em_names.extend(parts)

    ###########################
    ##Process band or wl ranges
    ###########################
    ##Check that we have band or wavelength windows
    # 检查两个参数 args.band_range 和 args.wl_range 是否同时被提供。如果这两个参数都包含了至少一个元素（即它们的长度都大于0），
    # 那么代码将抛出一个 RuntimeError 异常
    if (len(args.band_range) > 0) and (len(args.wl_range) > 0):
        raise RuntimeError("Can't mix wl_range and band_range entries")

    ##Convert to floats or ints
    # 处理程序的两个可选参数 args.band_range 和 args.wl_range，这两个参数分别代表不同的频率或波长范围。
    # 代码根据这两个参数的存在情况，分别解析它们的内容，并将解析后的结果存储在 band_ranges 或 wl_ranges 列表中
    if args.band_range:
        band_ranges = []
        wl_ranges = None
        for bstr in args.band_range:
            nobrackets = bstr.strip("[](){}")
            parts = nobrackets.split(",", maxsplit=2)
            if len(parts) < 2:
                raise RuntimeError(
                    "Could not convert band range to integers:" f" {bstr}"
                )
            band_ranges.append([int(p) for p in parts])

    elif args.wl_range:
        band_ranges = None
        wl_ranges = []
        for wstr in args.wl_range:
            nobrackets = wstr.strip("[](){}")
            parts = nobrackets.split(",", maxsplit=2)
            if len(parts) < 2:
                raise RuntimeError("Could not convert wl range to integers:" f" {wstr}")
            wl_ranges.append([float(p) for p in parts])
    else:
        raise RuntimeError(
            "Band or wavelength ranges must be specified with"
            "either --wl_range or --band_range"
        )

    if args.num_blocks:
        nblocks = []
        for nstr in args.num_blocks:
            nobrackets = nstr.strip("[](){}")
            parts = nobrackets.split(",", maxsplit=2)
            if len(parts) < 2:
                raise RuntimeError(
                    "Could not convert number of blocks to integers:" f" {nstr}"
                )
            nblocks.append([int(p) for p in parts])
    ##As a failsafe to prevent writing over one of the input files if command
    ## line argument order is mixed up, don't write output if the file exists
    if os.path.exists(args.output) and (not args.overwrite):
        raise RuntimeError(
            f"File {args.output} already exists. Check args"
            "with -v or use --overwrite if this is intentional"
        )

    #############################
    ##Create the AutoMCU instance
    #############################
    amcu = AutoMCU(
        em_csvs,
        em_counts,
        em_names=em_names,
        band_ranges=band_ranges,
        wl_ranges=wl_ranges,
        tied=args.tied,
        divide=args.divide,
        shade=args.shade,
        sum_to_one=args.sum_to_one,
        iterations=args.iterations,
        outscale=args.outscale,
        nodata=args.nodata,
        dtype=args.dtype,
        otype=args.otype,
        co=args.co,
        emfwhm=args.emfwhm,
        nointerp=args.nointerp,
        match_c=args.match_c,
    )

    ##########################
    ##Run AutoMCU on the image
    ##########################
    amcu.apply_image(
        args.input,
        args.scale,
        args.output,
        hdrpath=args.input_hdr,
        nblocks=args.num_blocks,
    )

    return


if __name__ == "__main__":
    main()
