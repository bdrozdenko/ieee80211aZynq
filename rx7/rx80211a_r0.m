function [srx, str, ffv, stx, stb, spi, spl, spq, ntf, cdi, cds] = rx80211a_r0()
[stx, stb, smf, spi, spl, spq, ntf, ~, ~, css] = tx80211a_r0();
% CONSTANTS
% Constant Number of data bits per frame
cdb = 24;
% Constant DeInterleaving index sequence
cdi = [0;3;6;9;12;15;18;21;24;27;30;33;36;39;42;45;1;4;7;10;13;16;19;22;25;28;31;34;37;40;43;46;2;5;8;11;14;17;20;23;26;29;32;35;38;41;44;47];
% Constant DeScrambling index sequence
cds = css;
% Constant Expected number of Frames
cef = uint16(ntf); % 36
% Constant number of samples per frame
cfs = 80;
% Constant number of transmitted samples
ctx = length(stx); % 3200
% FLAGS
% Flag: Preamble Found
fpf = false;
% Flag: Decode SIGNAL
fds = false;
% Flag: Decode MPDU
fdm = false;
% Flag: Last MPDU Frame Recovered
flf = false;
% INDICES (1-based)
% Index of Maximum Correlation
imc = 1;
% NUMBERS (counts)
% Number of decoded Bits (returned)
nb = 0;
% Number of preceding zero samples
np0 = randi([0,cfs-1],1,1);
% Number of received samples
nrs = ctx+80;
% Number of recovered MPDU frames
nrf = uint16(0);
% Number of Total recovered bits
ntb = 0;
% SEQUENCES (arrays)
% Sequence Receive Buffer: A FIFO for keeping received frames
srb = complex(zeros(3*cfs,1));
% Sequence Received: random # of zeros prepended to transmitted signal
srx = [complex(zeros(np0,1)); stx; complex(zeros(cfs-np0,1))];
% Sequence Received SIGNAL bits
srs = false(cdb,1);
% Sequence Received SERVICE+MPDU bits
srm = false(cdb*cef,1);
for ifr = 1:ceil(nrs/cfs)
    % Signal Received Frame (160 samples)
    srf(1:cfs,1) = srx((ifr-1)*cfs+1:(ifr*cfs),1);
    % Signal Receive Buffer (FIFO, 480 samples)
    srb(1:(2*cfs),1) = srb((cfs+1):(3*cfs),1); % Push back last 2 frames
    srb((2*cfs+1):(3*cfs)) = srf(1:cfs,1);   % Append latest received frame
    if ~fpf
        [fpf,fds,imc] = matchedfilter80211a(srb);
        if fpf
            fds = true;
            imc = imc+cfs;
        end
    elseif fds % SIGNAL  available, recover it now
        nrf = nrf + 1;  % Increment Rxd frames counter
        [nb,sb(1:3*cdb,1)] = recover80211aframe(srb(imc:(imc+cfs-1),1),nrf);
        fdm = true; fds = false;
    elseif fdm % SERVICE+MPDU available, recover both now
        nrf = nrf + 1;  % Increment Rxd frames counter
        [nb,sb(1:3*cdb,1)] = recover80211aframe(srb(imc:(imc+cfs-1),1),nrf);
        if (nb > 0)
            srm(ntb+1:ntb+nb,1) = sb(1:nb,1);
            ntb = ntb+nb;
        end
        if flf || (nrf>=cef)   % Upon recovery of last frame, cleanup
            fpf=false; fds=false; fdm=false; flf=false; 
            srb((cfs+1):(3*cfs),1) = complex(0.0); % Refresh Rx Buffer
        end
    else
        % Catch All: Program should never reach here
        error('Receiver entered an unexpected state. Exiting...');
    end
end
clear matchedfilter80211a recover80211aframe
if isequal(srm,stb)
    disp('PASS: Recovered Bits are equal to Transmitted Bits');
else
    disp('FAIL: Recovered Bits are NOT equal to Transmitted Bits');
end
[str,ffv] = unpack80211aframe(srm,smf);
end

function [fs1,fs0,imc] = matchedfilter80211a(srb)
% Matched Filter 80211a: Identifies presence/index of short&long preambles
% Inputs/Outputs:
% FPF: Flag Preamble Found, indicates PPDU Recovery Stage as enumeration: 
%   0: No preamble found
%   1: Short preamble found, but incomplete (<160 samples)
%   2: Short preamble found in full (160 samples)
%   3: Long preamble found, but incomplete (<160 samples)
%   4: Long preamble found in full (160 samples), decode SIGNAL next frame
%   5: Long preamble passed, SIGNAL ready to decode now (80 samples)
% IMC: Index of Maximum Correlation, identifies first valid sample
%   If FPF is 0     , IMC is 0
%   If FPF is 1 or 2, IMC is starting index of Short Preamble
%   If FPF is 3 or 4, IMC is starting index of Long Preamble
%   If FPF is 5     , IMC is starting index of SIGNAL
% SRF: Signal Receive Buffer: a FIFO of 480 complex samples from Receiver
if nargin < 1
    % Use transmitted sequence as received sequence
    tx80211a();
    srb = [complex(zeros(80,1)); stx; complex(zeros(80,1))];
end
persistent fp1 fp2 hfs hfl ixp
if isempty(fp2)
    fp2 = false;
    fp1 = false;
    ixp = 0;
    [sps, spl] = prepare80211apreamble();
    hfs = dsp.FIRFilter('Numerator',fliplr(sps'));
    hfl = dsp.FIRFilter('Numerator',fliplr(spl'));
end
% Constant number of Frame Samples
cfs = 80;
% Constant correlation Threshold for Definitive decision
ct3 = 2.0;
fp3 = false;
fs0 = false;
fs1 = false;
imc = 0;
% Actually run Matched Filter
sfs = step(hfs,srb);
sfl = step(hfl,srb);
[xcs,ixs] = max(sfs);
[xcl,ixl] = max(sfl);
xcs = abs(xcs); 
xcl = abs(xcl);
if fp1 || (xcl>ct3)
    fs1 = true; % Decode SIGNAL on Next (1) frame
    if fp1
        imc = ixp;
    else
        imc = mod(ixl,cfs)+1;
    end
    % Reset persistent variables
    fp2 = false; fp3 = false; ixp = 0;
elseif (xcs>ct3)
    fp3 = true; % Decode SIGNAL in Three frames
    ixp = mod(ixs,cfs)+1;
end
fp1 = fp2;
fp2 = fp3;
end

function [nb,sb] = recover80211aframe(srb,nfr)
% NFR: Frame Number: If 1, decode SIGNAL, if >1, decode SERVICE+MPDU
if nargin < 2
    nfr = 1;
end
persistent hbd hsc hvd ndf
% CONSTANTS
% Constant number of Coded Bits per OFDM symbol
ccb = uint16(48);
% Constant number of Cyclic Prefix samples
ccp = uint16(16);
% Constant number of Data Bits per OFDM symbol
cdb = 24;
% Constant Expected number of Frames
cef = uint16(36);
% Constant number of FFT bandwidths
cff = uint16(64);
% Constant Expected number of Bits
%ceb = cef*cdb;
% Constant Traceback Depth
ctd = 48;
% Constant number of Viterbi decoder delay Frames
cvf = ceil(ctd/cdb);
% FLAGS
% HANDLES TO SYSTEM OBJECTS
if isempty(hbd)
    hbd = comm.BPSKDemodulator('OutputDataType','double');
    hsc = comm.Scrambler('CalculationBase',2,'InitialConditions', ...
        [1,0,1,1,1,0,1],'Polynomial',[0,-4,-7]);
    hvd = comm.ViterbiDecoder('TrellisStructure',poly2trellis(7,[133,171]),...
        'InputFormat','Hard','TracebackDepth',ctd, ...
        'TerminationMethod','Continuous');
    ndf = 0;
end
% Indices for Bit Deinterleaving
ibi = zeros(48,1);
for k = 0:47, ibi(k+1)=(48/16)*mod(k,16)+floor(k/16); end
% Indices of Subcarriers (frequency bandwidths/channels)
is3 = [38:42,44:56,58:63,1:6,8:20,22:26].'+1; 
% Number of Bits (returned)
nb = 0;
% Sequence Bits (returned)
sb(1:(3*cdb),1) = false((3*cdb),1);
% Sequence Received OFDM demodulated
sro(1:cff,1) = fft(srb(ccp+1:ccp+cff,1),cff);
% Sequence Received BPSK deModulated
srm(1:ccb,1) = double(~step(hbd,sro(is3,1)));
% Sequence Received de-Interleaved
sri(1:ccb,1) = srm(ibi+1,1);
% Sequence Received Viterbi decoded
srv(1:cdb,1) = logical(step(hvd,sri));
ndf = ndf+1;
if (ndf > cvf)
    nb = cdb;
    if (ndf == cvf+1)
        sb(1:cdb,1) = srv(1:cdb,1);
    else
        % Scrambling sequence for seed 1011101, see Table G.15 p. 68
        sqs = step(hsc,zeros(24,1));
        % The Correct Descrambled Sequence
        sb(1:cdb,1) = xor(logical(sqs),srv(1:cdb,1)); % DESCRAMBLING
        if (nfr==cef)
            nb = nb+2*cdb;
            srv(1:cdb,1) = logical(step(hvd,zeros(ccb,1)));
            sqs = step(hsc,zeros(cdb,1));
            sb((cdb+1):(2*cdb),1) = xor(logical(sqs),srv(1:cdb,1));
            srv(1:cdb,1) = logical(step(hvd,zeros(ccb,1)));
            sqs = step(hsc,zeros(cdb,1));
            sb((2*cdb+1):(3*cdb),1) = xor(logical(sqs),srv(1:cdb,1));
            % Last 6 Bits are always 0s. Don't use Viterbi decoder output
            sb((3*cdb-5):(3*cdb),1) = logical(false(6,1));
            release(hsc); release(hvd); ndf = 0;
        end
    end
end
end

function [stx, stb, smf, spi, spl, spq, ntf, cii, cpp, css] = tx80211a_r0()
% Full Message Sequence
smf = sprintf(['Joy, bright spark of divinity,\nDaughter of Elysium,\n',...
    'Fire-insired we tread\nThy sanctuary.\nThy magic power re-unites',...
    '\nAll that custom has divided,\nAll men become brothers\n',...
    'Under the sway of thy gentle wings...\n']);
% PLCP PREAMBLE: Short & Long (each 160x1=8us)
[sps,spl] = prepare80211apreamble();
% SIGNAL (160x1=4us) + SERVICE+MPDU
[ntf,sts,stb] = prepare80211aframe(smf);
% Full Transmit Sequence
sp  = [sps; spl];
spf = fi(sp,true,16,14,'RoundingMethod','Round','OverflowAction','Saturate');
spi = real(spf);
spq = imag(spf);
stx = [sp; sts(:)];
%% Constants
% Constant Interleaving Indices
cii = uint8(reshape(reshape(0:47,16,3).',48,1));
% Constant Pilot Polarities
cpp = [1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,1,1,1,1,1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,1,-1,-1,1,-1,-1,1,1,1,1,1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,-1,-1,-1,-1,-1].';
% Constant Scrambling Sequence
css = repmat(logical([0;1;1;0;1;1;0;0;0;0;0;1;1;0;0;1;1;0;1;0;1;0;0;1;1;1;0;0;1;1;1;1;0;1;1;0;1;0;0;0;0;1;0;1;0;1;0;1;1;1;1;1;0;1;0;0;1;0;1;0;0;0;1;1;0;1;1;1;0;0;0;1;1;1;1;1;1;1;0;0;0;0;1;1;1;0;1;1;1;1;0;0;1;0;1;1;0;0;1;0;0;1;0;0;0;0;0;0;1;0;0;0;1;0;0;1;1;0;0;0;1;0;1;1;1;0;1]),2,1);
end

function [ssp, slp] = prepare80211apreamble()
% Constant number of Preamble (either short or long) Samples
cps = uint16(160);
% Constant number of FFT/IFFT points: always 64
cff = uint16(64);
% Constant number of Used subCarriers: Constant 52 + 1 for DC
cuc = uint16(53);
% Number of Long Training Symbols: Constant 2
ntl = uint16(2);
% Number of Short Training Symbols: Constant 10
nts = uint16(10);
% Preallocate memory for vectors
% Sequence Short preamble Training symbols
sst = complex(zeros(cuc,1)); 
% Sequence Short preamble OFDM modulated
sso = complex(zeros(cff,1));
% Sequence Short Preamble (truncated to 160 samples)
ssp = complex(zeros(cps,1));
% Sequence Long Preamble Training symbols
slt = complex(zeros(cuc,1)); 
% Sequence Long Preamble OFDM modulated
slo = complex(zeros(cff,1));
% Sequence Long Preamble complete
slp = complex(zeros(cff*ntl+cff/2,1));
% The 3 components that comprise a short training symbol
st0 = [0;0;0;0]; st1 = [0;0;1+1i;0]; st2 = [0;0;-1-1i;0]; 
% Short Training Symbol: 12 subcarriers used between -26 and 26
sst(1:cuc,1) = sqrt(13/6)*[st1;st2;st1;st2;st2;st1; st0; ...
    st2;st2;st1;st1;st1;st1;0];
% Short Training Symbol, OFDM Modulated
sso(1:cff,1) = ifft([sst(27:53,1);complex(zeros(11,1));sst(1:26,1)],64);
% Sequence Short Preamble: 10 Short OFDM-Modulated Training Symbols
ssp(1:cps,1) = [sso; sso; sso(1:cff/2)];

% Long Training Symbol: 52 subcarriers plus a zero value at DC, 
%                    modulated by the elements of the sequence L, given by: 
slt(1:cuc,1) = complex([1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;-1;-1;1;1;-1;1;-1;1;1;1;1; ...
    0; 1;-1;-1;1;1;-1;1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;1;-1;1;1;1;1], ...
    zeros(cuc,1));
% Long Training Symbol, OFDM Modulated
slo(1:cff,1) = ifft([slt(27:53,1);complex(zeros(11,1));slt(1:26,1)],64);
% Sequence Preamble Long: 2 Long OFDM-Modulated Training Symbols
slp((cff/2+1):(cff*ntl+cff/2),1) = repmat(slo,ntl,1);
% Prepend OFDM-Modulated Long Training Symbols w/ 16-symbol Cyclic Prefix (GI2)
slp(1:(cff/2),1) = slo(end-31:end,1);
end

function [ntf,sts,stb] = prepare80211aframe(smf)
if nargin < 1
    % Sequence Message in Full
    smf = sprintf(['Joy, bright spark of divinity,\nDaughter of Elysium,\n',...
        'Fire-insired we tread\nThy sanctuary.\nThy magic power re-unites',...
        '\nAll that custom has divided,\nAll men become brothers\n',...
        'Under the sway of thy gentle wings...\n']);
end
% CONSTANTS
% Bit Ordering
cml = 'right-msb'; 
% Constant number of Coded bits per OFDM symbol (see Table 78 p. 9)
cco = 48;
% Constant rate 1/2 Convolutional coder Upsampling factor
c12 = 1/2;
% Constant rate 3/4 Convolutional coder Upsampling factor
c34 = 3/4;
% Constant number of SIGNAL data Bits per OFDM symbol (see Table 78 p. 9)
csb = cco*c12;
% Constant number of SVC+MPDU data Bits per OFDM symbol (see Table 78 p. 9)
cmb = cco*c34;
% Constant number of Coded bits per Subcarrier/BPSK symbol (Table 78 p. 9)
%ccs = 1;
% Constant Polarity of Pilot subcarriers
cpp = [1,1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,1,1,1,1,1,1,-1,1,...
    1,1,-1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,1,-1,-1,1,-1,-1,1,1,1,1,1,-1,-1,1,1,...
    -1,-1,1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,1,-1,1,...
    -1,-1,-1,-1,-1,1,-1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,-1,1,1,1,-1,-1,-1,-1,-1,-1,-1];
% HANDLES TO SYSTEM OBJECTS
hbm = comm.BPSKModulator('OutputDataType','double');
hcg = comm.CRCGenerator([32,26,23,22,16,12,11,10,8,7,5,4,2,1,0],...
        'InitialConditions',1,'DirectMethod',true,...
        'ReflectInputBytes',true,'ReflectChecksums',false,'FinalXOR',true);
hif = dsp.IFFT('Normalize',true,...
    'FFTLengthSource','Property','FFTLength',64);
hsc = comm.Scrambler('CalculationBase',2,'InitialConditions', ...
    [1,0,1,1,1,0,1],'Polynomial',[0,-4,-7]);
h12 = comm.ConvolutionalEncoder('TrellisStructure',poly2trellis(7,[133,171]),...
    'TerminationMethod','Continuous');
h34 = comm.ConvolutionalEncoder('TrellisStructure',poly2trellis(7,[133,171]),...
    'TerminationMethod','Continuous','PuncturePatternSource','Property',...
    'PuncturePattern',[1;1;1;0;0;1]);
% INDICES
% Indices for Bit Interleaving
ibi = zeros(48,1);
for k = 0:47, ibi(k+1)=(48/16)*mod(k,16)+floor(k/16); end
% Indices of Subcarriers & Pilots (frequency bandwidths/channels)
is0 = [0:5,32,59:63].'+1;                       % Indices of Unused Channels
is1 = [6:10,12:24,26:31,33:38,40:52,54:58].'+1; % Indices of Coded Bit Subcarriers
is2 = [11,25,39,53].'+1;                        % Indices of Pilot Subcarriers
% INTEGER NUMBERS
% Number of Cyclic Prefix samples
ncp = uint16(16);
% Number of FFT bandwidths
nff = uint16(64);
% Number of Transmitted Samples per frame
nts = ncp+nff;
% Number of Header Characters
nhc = 24;
% Number of Header Bits
nhb = nhc*8;
% Number of Message Characters
nmc = 72;
% Number of Message Bits
nmb = nmc*8;
% Number of CRC (FCS) Characters
ncc = 4;
% Number of CRC (FCS) Bits
ncb = ncc*8;
% Number of DATA (total MPDU) Characters
%ndc = nhc+nmc+ncc;
% Number of DATA (total MPDU) Bits
ndb = nhb+nmb+ncb;
% Number of serVice Bits
nvb = 16;
% Number of taiL Bits
nlb = 6;
% Number of service+DATA+tail Bits
nxb = ndb+nvb+nlb;
% Number of service+DATA+tail+pad OFDM symbols / frames
nzo = ceil(nxb/csb);
% Number of service+DATA+tail+pad Bits
nzb = nzo*csb;
% Number of Pad bits
npb = nzb - nxb;
% Number of Transmission Frames
ntf = nzo+1;

% SEQUENCES (binary, hex, character, etc.)
% PLCP HEADER: SIGNAL field
ssr = logical([1;0;1;1]);                   % RATE
ssl = logical([0;0;1;0;0;1;1;0;0;0;0;0]);   % LENGTH
sst = logical([0;0;0;0;0;0]);               % TAIL
ssb = logical([ssr; false; ssl; false; sst]);
% MAC Header Sequence (in hex)
shh = '0402002e006008cd37a60020d6013cf1006008ad3baf0000';   %HEXADECIMAL
shd = hex2dec(reshape(shh,2,nhc).');                        %DECIMAL
shb(1:nhb,1) = reshape(de2bi(shd,8,cml).',nhb,1);           %BINARY
% Truncated Message Sequence
smc = smf(1:nmc);                                           %CHARACTERS
smd = double(smc);                                          %DECIMAL
%smh = dec2hex(smd);                                         %HEXADECIMAL
smb(1:nmb,1) = reshape(de2bi(smd,8,cml).',nmb,1);           %BINARY
% PSDU/MPDU: MAC Header + Message Body
sda = [shb;smb];
% Message to Transmit (with Appended CRC)
sdb = step(hcg,sda);                    % BINARY
% Known problem: CRC is left-msb--reverse bit order for each byte
sdb(end-31:end,1) = [sdb(end-24:-1:end-31,1);sdb(end-16:-1:end-23,1); ...
    sdb(end-8:-1:end-15,1);sdb(end:-1:end-7,1)];
%sdd = bi2de(reshape(sdb,8,ndc).',cml);  % DECIMAL
%sdc = char(sdd.');                      % CHARACTERS
%sdh = dec2hex(sdd);                     % HEXADECIMAL
% Sequence of SERVICE bits
svb = false(nvb,1);
% Sequence of TAIL bits
slb = false(nlb,1);
% Sequence of PAD bits
spb = false(npb,1);
% Sequence of SERVICE+MPDU+TAIL+PAD bits
szb = [svb;logical(sdb);slb;spb];
% Sequence of Transmission bits
stb = [ssb;szb];
% PREALLOCATION OF MEMORY FOR ARRAYS
stc = false(csb,ntf);
ste = zeros(cco,ntf);
sti = zeros(cco,ntf); 
stm = complex(zeros(nff,ntf));
sto = complex(zeros(nts,ntf));
sts = complex(zeros(nts,ntf));

% PREPARE SERVICE+MPDU
for itf = 1:ntf
    if (itf == 1)   % SIGNAL: No Scrambling
        stc(1:csb,itf) = stb((itf-1)*csb+1:(itf*csb),1);
        % Sequence of Transmission convolutionally Encoded bits
        ste(1:cco,itf) = step(h12,stc(1:csb,itf));          % RATE 1/2 CONVOLUTIONAL ENCODING
    else            % SERVICE+MPDU: Scrambling
        % Scrambling sequence for seed 1011101, see Table G.15 p. 68
        sqs = step(hsc,zeros(24,1));
        % The Correct Scrambling Sequence, Table G.16 p. 69
        stc(1:csb,itf) = xor(logical(sqs),stb((itf-1)*csb+1:(itf*csb),1));
        % Sequence of Transmission convolutionally Encoded bits
        ste(1:cco,itf) = step(h12,stc(1:csb,itf));          % RATE 1/2 CONVOLUTIONAL ENCODING
    end
    % Sequence of Transmission Interleaved coded bits
    sti(ibi+1,itf) = ste(1:cco,itf);                    % BLOCK INTERLEAVING
    % Sequence of Transmission BPSK Modulated and mapped symbols
    stm(is0,itf) = complex(0.0);
    stm(is1,itf) = complex(step(hbm,~sti(1:cco,itf)));
    stm(is2,itf) = complex(cpp(itf)*[1.0;1.0;1.0;-1.0]);% PILOT SUBCARRIERS
    % Sequence of OFDM modulated symbols
    sto(ncp+1:nts,itf) = step(hif,[stm(33:64,itf);stm(1:32,itf)]);
    %ifft([stm(33:64,itf);stm(1:32,itf)],nff);
    sto(1:ncp,itf) = sto(nff+1:nts,itf);
    % Sequence of Transmission (Time) Samples
    sts(1:nts,itf) = sto(1:nts,itf);
end
% CLEANUP: Release System objects
release(hcg); release(hsc); release(h12);
save txol stb stc ste sti stm sto sts
end

function [str,ffv] = unpack80211aframe(sb,smf)
if nargin < 1
    % Sequence Message in Full
    smf = sprintf(['Joy, bright spark of divinity,\nDaughter of Elysium,\n',...
        'Fire-insired we tread\nThy sanctuary.\nThy magic power re-unites',...
        '\nAll that custom has divided,\nAll men become brothers\n',...
        'Under the sway of thy gentle wings...\n']);
end
% HANDLES to System objects
hcd = comm.CRCDetector([32,26,23,22,16,12,11,10,8,7,5,4,2,1,0],...
    'InitialConditions',1,'DirectMethod',true,...
    'ReflectInputBytes',true,'ReflectChecksums',false,'FinalXOR',true);
str.aSIGNAL.a       = sb(1:24,1);
str.aSIGNAL.bRATE   = sb(1:4,1);
str.aSIGNAL.cLENGTH = sb(6:17,1);
str.aSIGNAL.dTAIL   = sb(19:24,1);
str.bSERVICE        = sb(24+(1:16),1);
str.cMPDU.a         = sb(40+(1:800),1);
str.cMPDU.bMAChdr.a = sb(40+(1:192),1);
str.cMPDU.bMAChdr.bFrameCtrl.a                  = sb(40+(1:16),1);
str.cMPDU.bMAChdr.bFrameCtrl.bProtocolVersion   = sb(40+(1:2),1);
str.cMPDU.bMAChdr.bFrameCtrl.cType              = sb(42+(1:2),1);
str.cMPDU.bMAChdr.bFrameCtrl.dSubtype           = sb(44+(1:4),1);
str.cMPDU.bMAChdr.bFrameCtrl.eToDS              = sb(49,1);
str.cMPDU.bMAChdr.bFrameCtrl.fFromDS            = sb(50,1);
str.cMPDU.bMAChdr.bFrameCtrl.gMoreFrag          = sb(51,1);
str.cMPDU.bMAChdr.bFrameCtrl.hRetry             = sb(52,1);
str.cMPDU.bMAChdr.bFrameCtrl.iPwrMgmt           = sb(53,1);
str.cMPDU.bMAChdr.bFrameCtrl.jMoreData          = sb(54,1);
str.cMPDU.bMAChdr.bFrameCtrl.kWEP               = sb(55,1);
str.cMPDU.bMAChdr.bFrameCtrl.lOrder             = sb(56,1);
str.cMPDU.bMAChdr.cDurationID.a                 = sb(56+(1:16),1);
str.cMPDU.bMAChdr.cDurationID.bValue            = bi2de(sb(56+(1:16),1));
str.cMPDU.bMAChdr.dAddress1                     = sb(72+(1:48),1);
str.cMPDU.bMAChdr.eAddress2                     = sb(120+(1:48),1);
str.cMPDU.bMAChdr.fAddress3                     = sb(168+(1:48),1);
str.cMPDU.bMAChdr.gSequenceCtrl.a               = sb(216+(1:16),1);
str.cMPDU.bMAChdr.gSequenceCtrl.bFragmentNum    = sb(216+(1:4),1);
str.cMPDU.bMAChdr.gSequenceCtrl.cSequenceNum    = sb(220+(1:12),1);
str.cMPDU.bMAChdr.hHex = dec2hex(bi2de(reshape(sb(40+(1:192),1),8,24).'),2);
str.cMPDU.cFrameBody.a      = sb(232+(1:576),1);
str.cMPDU.cFrameBody.bDec   = bi2de(reshape(sb(232+(1:576),1),8,72).');
str.cMPDU.cFrameBody.cChar  = char(str.cMPDU.cFrameBody.bDec.');
str.cMPDU.cFrameBody.dHex   = dec2hex(str.cMPDU.cFrameBody.bDec,2);
str.cMPDU.dFCS.a    = sb(808+(1:32),1);
str.cMPDU.dFCS.bHex = dec2hex(bi2de(reshape(sb(808+(1:32),1),8,4).'),2);
str.dTAIL           = sb(840+(1:6),1);
str.ePAD            = sb(846+(1:18),1);
% Check CRC
smp = str.cMPDU.a;
% Known problem: CRC is left-msb--reverse bit order for each byte
smp(end-31:end,1) = [smp(end-24:-1:end-31,1);smp(end-16:-1:end-23,1); ...
    smp(end-8:-1:end-15,1);smp(end:-1:end-7,1)];
[~,fer] = step(hcd,smp);
if fer
    disp('FAIL: Frame Check Sequence (FCS) yielded a CRC error.');
else
    disp('PASS: Frame Check Sequence (FCS) yielded the correct CRC.');    
end
% Check Message Characters
if isequal(str.cMPDU.cFrameBody.cChar,smf(1:72))
    disp('PASS: Recovered message is equal to Transmitted message.');
else
    disp('FAIL: Recovered message is NOT equal to Transmitted message.');
end
ffv = ~fer;
end
