import bisect
import collections
import errno
import os
import sys
import timeit
from enum import Enum, EnumMeta, unique
from io import StringIO

import natsort as ns
import numpy as np
import pandas as pd
import pyranges as pr
from rich.console import Console

# write console messages to standard error
console = Console(stderr=True)

APP_NAME = "filter-regions"
MAX_ELEMS = 100_000


class FilterMeta(EnumMeta):
    def __contains__(cls, item):
        return item in [v.value for v in cls.__members__.values()]


@unique
class FilterMethod(str, Enum, metaclass=FilterMeta):
    pq = "pq"
    wis = "wis"
    maxmean = "maxmean"


class FilterMethodDescription(str, Enum, metaclass=FilterMeta):
    pq = "Priority-Queue (PQ)"
    wis = "Weighted Interval Scheduler (WIS)"
    maxmean = "Max-Mean Sweep (MaxMean)"


@unique
class FilterInputType(str, Enum, metaclass=FilterMeta):
    vector = "vector"
    bedgraph = "bedgraph"


@unique
class FilterAggregationMethod(str, Enum, metaclass=FilterMeta):
    min = "min"
    max = "max"
    mean = "mean"
    sum = "sum"
    median = "median"
    variance = "variance"
    percentile = "percentile"


class Filter:
    def __init__(
        self,
        method,
        input,
        input_type,
        window_bins=None,
        bin_size=200,
        exclusion_size=24800,
        preserve_cols=False,
        aggregation_method=FilterAggregationMethod.max,
        percentile=0.95,
        max_elements=MAX_ELEMS,
        quiet=False,
    ) -> None:
        self.method: str = method
        self.input: str = input
        self.input_type: str = input_type
        self.window_bins: int = window_bins if window_bins else None
        self.bin_size: int = bin_size
        self.exclusion_size: int = exclusion_size
        self.input_df: pd.DataFrame = None
        self.output_df: pd.DataFrame = None
        self.preserve_cols: bool = preserve_cols
        self.aggregation_method: str = aggregation_method
        self.percentile: float = percentile
        self.max_elements: int = max_elements
        self.quiet: bool = quiet
        # ignore SettingWithCopyWarning
        # cf. https://towardsdatascience.com/explaining-the-settingwithcopywarning-in-pandas-ebc19d799d25
        pd.options.mode.chained_assignment = None

    def input_df(self, df: pd.DataFrame) -> None:
        self.input_df = df

    def output_df(self, df: pd.DataFrame) -> None:
        self.output_df = df

    def preserve_cols(self, flag: bool) -> None:
        self.preserve_cols = flag

    def input_type(self, it: str) -> None:
        self.input_type = it

    def max_elements(self, me: int) -> None:
        self.max_elements = me

    def read(self) -> None:
        df = None
        if not self.quiet:
            console.print(
                f"[bold blue]{APP_NAME}[/] | [bold green]Run[/] | Reading input..."
            )
        if isinstance(self.input, str):
            if not os.path.exists(self.input):
                console.print(
                    f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Must specify valid input path \[{self.input}]"
                )
                sys.exit(errno.EINVAL)
            if self.method == "pq":
                df = pd.read_csv(
                    self.input,
                    sep="\t",
                    header=None
                )
            elif self.method == "maxmean":
                df = pd.read_csv(
                    self.input,
                    sep="\t",
                    header=None
                )
            elif self.method == "wis":
                if self.input_type == "bedgraph":
                    df = pr.read_bed(self.input)
                elif self.input_type == "vector":
                    df = pd.read_csv(
                        self.input,
                        sep="\t",
                        header=None
                    )
        elif isinstance(self.input, np.ndarray):
            s = self.input.shape
            l = len(s)
            if l != 1 and l != 2:
                console.print(
                    f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Numpy array is not a 1D vector or 2D matrix \[{l} != 1 or 2]"
                )
                sys.exit(errno.EINVAL)
            if l == 2 and s[1] != 4:
                console.print(
                    f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Numpy matrix has incorrect number of columns \[{s[1]} != 4]"
                )
                sys.exit(errno.EINVAL)
            if self.method == "pq":
                df = pd.DataFrame(self.input)
                if l == 1:
                    df = df.transpose()
            elif self.method == "maxmean":
                df = pd.DataFrame(self.input)
                if l == 1:
                    df = df.transpose()
            elif self.method == "wis":
                df = pd.DataFrame(self.input)
                if l == 1:
                    df = df.transpose()
        if self.input_type == "bedgraph":
            df.columns = ["Chromosome", "Start", "End", "Score"]
        elif self.input_type == "vector":
            df = df.rename({0: "Score"}, axis='index')
            df = df.transpose().reset_index(drop=True)
            df["Start"] = df.index
            df["End"] = df["Start"] + 1
            if self.method == "wis":
                df["Chromosome"] = "chrSentinel"
                df = df[["Chromosome", "Start", "End", "Score"]]
            else:
                df = df[["Start", "End", "Score"]]
        df["Score"] = df["Score"].astype(float) # req'd for sort to work correctly
        if self.method == "wis" and self.input_type == "vector":
            df = pr.PyRanges(df)
        self.input_df = df
        self.window_bins = (
            self.window_bins
            if self.window_bins
            else (self.bin_size + self.exclusion_size) // self.bin_size
        )  # window size (in units of bins, not nt)

    def write(self, output: None) -> None:
        o = StringIO()
        try:
            self.output_df.to_csv(o, sep="\t", index=False, header=False)
            if not output:
                sys.stdout.write(f"{o.getvalue()}")
            else:
                with open(output, "wb") as ofh:
                    ofh.write(f"{o.getvalue()}".encode("utf-8"))
        except AttributeError as err:
            console.print(
                f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Could not write output \[{err}]"
            )
            sys.exit(errno.EINVAL)

    def output_df(self) -> pd.DataFrame:
        return self.output_df

    def filter(self) -> None:
        try:
            if not self.quiet:
                console.print(
                    f"[bold blue]{APP_NAME}[/] | [bold green]Run[/] | Filtering with {self.method} method..."
                )
                console.print(
                    f"[bold blue]{APP_NAME}[/] | [bold green]Run[/] | Aggregating score window with {self.aggregation_method} method..."
                )
            start = timeit.default_timer()
            run = getattr(self, self.method)
            df = run()
            if self.aggregation_method == "max":
                df.Score = df.RollingMax
            elif self.aggregation_method == "min":
                if self.preserve_cols:
                    df.Score = df.RollingMin
                else:
                    console.print(
                        f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Cannot access aggregation score without --preserve-cols"
                    )
                    sys.exit(errno.EINVAL)
            elif self.aggregation_method == "mean":
                df.Score = df.RollingMean
            elif self.aggregation_method == "sum":
                if self.preserve_cols:
                    df.Score = df.RollingSum
                else:
                    console.print(
                        f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Cannot access aggregation score without --preserve-cols"
                    )
                    sys.exit(errno.EINVAL)
            elif self.aggregation_method == "median":
                if self.preserve_cols:
                    df.Score = df.RollingMedian
                else:
                    console.print(
                        f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Cannot access aggregation score without --preserve-cols"
                    )
                    sys.exit(errno.EINVAL)
            elif self.aggregation_method == "variance":
                if self.preserve_cols:
                    df.Score = df.RollingVariance
                else:
                    console.print(
                        f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Cannot access aggregation score without --preserve-cols"
                    )
                    sys.exit(errno.EINVAL)
            elif self.aggregation_method == "percentile":
                if self.preserve_cols:
                    df.Score = df.RollingPercentile
                else:
                    console.print(
                        f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Cannot access aggregation score without --preserve-cols"
                    )
                    sys.exit(errno.EINVAL)
            self.output_df = df
            end = timeit.default_timer()
            if not self.quiet:
                console.print(
                    f"[bold blue]{APP_NAME}[/] | [bold green]Run[/] | Method completed in: {(end - start):.2f} s"
                )
                console.print(
                    f"[bold blue]{APP_NAME}[/] | [bold green]Run[/] | Method found: {len(self.output_df.index)} of {len(self.input_df.index)} ({(len(self.output_df.index) / len(self.input_df.index)):.4f}) elements"
                )
        except TypeError:
            console.print(
                f"[bold blue]{APP_NAME}[/] | [bold red]Error[/] | Method could not be applied"
            )
            sys.exit(errno.EINVAL)

    def pq(self) -> pd.DataFrame:
        return self.maxmean(pq=True)

    def wis(self) -> pd.DataFrame:
        d = self.input_df
        m = (
            pr.PyRanges(d).merge()
            if isinstance(d, pd.DataFrame)
            else d.merge()
        )
        df = (
            d
            if isinstance(d, pd.DataFrame)
            else d.as_df()
        )
        self.input_df = df
        df = df.rename(columns={"Name": "Score"})
        # update the start and end loci based on window size (bin units)
        df["Start"] = df.Start.shift(self.window_bins // 2)
        df["End"] = (
            df.End.shift(-(self.window_bins // 2))
            if self.window_bins % 2
            else df.End.shift(-(self.window_bins // 2 - 1))
        )
        df = df.dropna()
        df = df[df["Start"] < df["End"]]
        df["Start"] = df["Start"].astype(int)
        df["End"] = df["End"].astype(int)
        # intervals in the df must be sorted for the trace to work correctly (natural sort on chromosome)
        df["OriginalIdx"] = df.index
        df = df.sort_values(by=["Chromosome", "Start", "End"])
        df = df.reindex(
            index=ns.order_by_index(
                df.index, ns.index_natsorted(zip(df.Chromosome, df.Start, df.End))
            )
        )
        df = df.reset_index(drop=True)
        # get rolling summary statistics
        df["RollingMin"] = df["Score"].rolling(self.window_bins, center=True).min()
        df["RollingMax"] = df["Score"].rolling(self.window_bins, center=True).max()
        df["RollingMean"] = df["Score"].rolling(self.window_bins, center=True).mean()
        df["RollingSum"] = df["Score"].rolling(self.window_bins, center=True).sum()
        df["RollingMedian"] = df["Score"].rolling(self.window_bins, center=True).median()
        df["RollingVariance"] = df["Score"].rolling(self.window_bins, center=True).var()
        df["RollingPercentile"] = df["Score"].rolling(self.window_bins, center=True).quantile(self.percentile)
        # get rid of edges
        df = df.dropna().reset_index(drop=True)
        df["MethodIdx"] = df.index
        n = len(df.index)
        # build accumulator-based dict for calculation of absolute coordinates
        acc_chr = m.as_df().loc[:, "Chromosome"].copy().values
        acc_abs = np.concatenate([[0], m.as_df().loc[:, "End"].copy().values[:-1]])
        j = len(acc_abs) - 1
        while j > 0:
            acc_abs[j] = np.sum(acc_abs[0 : j + 1])
            j -= 1
        acc = dict(zip(acc_chr, acc_abs))
        # for every interval j, compute the rightmost disjoint interval i, where i < j.
        s = df.loc[:, "Start"].copy().values
        e = df.loc[:, "End"].copy().values
        # translate per-chr coord space to abs coord space
        c = df.loc[:, "Chromosome"].copy().values
        translate_genomic_coords_to_abs = np.vectorize(lambda x, y: y + acc[x])
        s = translate_genomic_coords_to_abs(c, s)
        e = translate_genomic_coords_to_abs(c, e)
        # look for nearest disjoint interval to left of current interval
        p = []
        for j in range(n):
            i = bisect.bisect_right(e, s[j]) - 1
            p.append(i)
        # set up initial opt(imum) table values.
        opt = collections.defaultdict(int)
        opt[-1] = 0
        opt[0] = 0
        # forwards pass to get best-scoring path
        for j in range(1, n):
            score = df.loc[df.index[j], "Score"]
            opt[j] = max(score + opt[p[j]], opt[j - 1])

        # backwards trace to retrieve path
        q = []
        j = n - 1
        while j >= 0:
            score = df.loc[df.index[j], "Score"]
            if score + opt[p[j]] > opt[j - 1]:
                q.append(j)
                j = p[j]  # jump to the nearest disjoint interval to the left
            else:
                j -= 1  # try the "next" interval, one to the left
        # sort indices of qualifiers
        q.sort()
        df = df.iloc[q]
        if self.input_type == "vector":
            df = df[["Start", "End", "Score"]] if not self.preserve_cols else df.iloc[: , 1:]
        if len(df.index) > self.max_elements:
            df = df \
                .sort_values(by=["Score"], ascending=False) \
                .head(self.max_elements) \
                .sort_index(axis=0)
        return df

    def maxmean(self, pq: bool = False) -> pd.DataFrame:
        df = self.input_df
        # save original indices
        df["OriginalIdx"] = df.index
        # update the start and end loci based on window size (bin units)
        df["Start"] = df.Start.shift(self.window_bins // 2)
        df["End"] = (
            df.End.shift(-(self.window_bins // 2))
            if self.window_bins % 2
            else df.End.shift(-(self.window_bins // 2 - 1))
        )
        df = df.dropna()
        df["Start"] = df["Start"].astype(int)
        df["End"] = df["End"].astype(int)
        df = df.reset_index(drop=True)
        # get rolling summary statistics
        df["RollingMax"] = df["Score"].rolling(self.window_bins, center=True).max()
        df["RollingMean"] = df["Score"].rolling(self.window_bins, center=True).mean()
        if self.preserve_cols:
            df["RollingMin"] = df["Score"].rolling(self.window_bins, center=True).min()
            df["RollingSum"] = df["Score"].rolling(self.window_bins, center=True).sum()
            df["RollingMedian"] = df["Score"].rolling(self.window_bins, center=True).median()
            df["RollingVariance"] = df["Score"].rolling(self.window_bins, center=True).var()
            df["RollingPercentile"] = df["Score"].rolling(self.window_bins, center=True).quantile(self.percentile)
        # get rid of edges
        df = df.dropna().reset_index(drop=True)
        # drop regions which fall over two chromosomes
        if self.input_type == "bedgraph":
            df = df.drop(
                np.where(df.iloc[:, 1].to_numpy() >= df.iloc[:, 2].to_numpy())[0]
            ).reset_index(drop=True)
        # save new indices before sorting
        df["MethodIdx"] = df.index
        # sort by max, by mean, and then by score, unless we are
        # using the pq method, in which case we simply prioritize
        # on the original score column
        df = (
            df.sort_values(by=["RollingMax", "RollingMean", "Score"], ascending=False)
            if not pq
            else df.sort_values(by=["Score"], ascending=False)
        )
        # run filter algorithm
        n = len(df)
        hits = np.zeros(n, dtype=bool)
        indices = []
        k = self.max_elements
        # precalculate values for regions from individual indices
        df["MethodStartIdx"] = df["MethodIdx"] - self.window_bins // 2
        df["MethodEndIdx"] = df["MethodIdx"] + self.window_bins // 2 + 1 if self.window_bins % 2 else df["MethodIdx"] + self.window_bins // 2
        df.loc[df["MethodStartIdx"] < 0, "MethodStartIdx"] = 0
        df.loc[df["MethodEndIdx"] > n, "MethodEndIdx"] = n
        loc_MethodIdx = df.columns.get_loc("MethodIdx")
        loc_MethodStartIdx = df.columns.get_loc("MethodStartIdx")
        loc_MethodEndIdx = df.columns.get_loc("MethodEndIdx")
        for loc in range(n):
            if k <= 0:
                break
            start = df.iloc[loc, loc_MethodStartIdx]
            stop = df.iloc[loc, loc_MethodEndIdx]
            if not np.any(hits[start:stop]):
                hits[start:stop] = True
                indices.append(df.iloc[loc, loc_MethodIdx])
                k -= 1
        df = df.loc[indices]
        # resort input by original order
        df = df.sort_values(by=["OriginalIdx"], ascending=True)
        # clip columns, unless we preserve them
        if not self.preserve_cols:
            if self.input_type == "bedgraph":
                df = df[["Chromosome", "Start", "End", "Score"]]
            elif self.input_type == "vector":
                df = df[["Start", "End", "Score"]]
        df = df.reset_index(drop=True)
        return df
