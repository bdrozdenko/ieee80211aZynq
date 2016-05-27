function [stx, stb, spi, spq, ntf, cii, cpp, css] = tx80211a_r0()
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
