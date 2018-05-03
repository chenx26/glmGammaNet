function (data, ..., max.freq = 0.5, twosided = FALSE, keep = 1)
{
        data.fft = fft(data)
        N = length(data)
        tmp = Mod(data.fft[2:floor(N/2)])^2/N
        tmp = sapply(tmp, function(x) max(1e-05, x))
        freq = ((1:(floor(N/2) - 1))/N)
        tmp = tmp[1:floor(length(tmp) * keep)]
        freq = freq[1:floor(length(freq) * keep)]
        if (twosided) {
                tmp = c(rev(tmp), tmp)
                freq = c(-rev(freq), freq)
        }
        return(list(spec = tmp, freq = freq))
}
