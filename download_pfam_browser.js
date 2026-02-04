// =============================================================================
// HMMER/PFAM Bulk Results Downloader - Browser Console Script
// =============================================================================
//
// INSTRUCTIONS:
// 1. Open your HMMER results page in Chrome/Firefox:
//    https://www.ebi.ac.uk/Tools/hmmer/results/50b4593e-293e-4875-b211-1454cba3cbb8/score
// 2. Open Developer Tools (F12 or Ctrl+Shift+I)
// 3. Go to the "Console" tab
// 4. Copy and paste this entire script
// 5. Press Enter to run
//
// The script will download all results as JSON files.
// =============================================================================

(async function downloadAllPfamResults() {
    const JOB_ID = '50b4593e-293e-4875-b211-1454cba3cbb8';
    const BASE_API = 'https://www.ebi.ac.uk/Tools/hmmer/api/v1';

    console.log('Starting HMMER results download...');
    console.log('Job ID:', JOB_ID);

    // Helper function to download JSON as file
    function downloadJSON(data, filename) {
        const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    // Helper function to download text as file
    function downloadText(text, filename) {
        const blob = new Blob([text], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    // Helper function to delay between requests
    function delay(ms) {
        return new Promise(resolve => setTimeout(resolve, ms));
    }

    try {
        // 1. Fetch main results
        console.log('Fetching main results...');
        const resultsResponse = await fetch(`${BASE_API}/result/${JOB_ID}`, {
            headers: { 'Accept': 'application/json' }
        });

        if (!resultsResponse.ok) {
            throw new Error(`Failed to fetch results: ${resultsResponse.status}`);
        }

        const results = await resultsResponse.json();
        console.log('Main results fetched successfully');

        // Download main results
        downloadJSON(results, 'hmmer_all_results.json');
        console.log('Downloaded: hmmer_all_results.json');

        await delay(500);

        // 2. Fetch domain architecture
        console.log('Fetching domain architecture...');
        try {
            const archResponse = await fetch(`${BASE_API}/architecture/${JOB_ID}`, {
                headers: { 'Accept': 'application/json' }
            });

            if (archResponse.ok) {
                const archData = await archResponse.json();
                downloadJSON(archData, 'hmmer_domain_architectures.json');
                console.log('Downloaded: hmmer_domain_architectures.json');
            }
        } catch (e) {
            console.log('Architecture data not available:', e.message);
        }

        await delay(500);

        // 3. Extract and process hits
        let hits = [];
        if (results.results && results.results.hits) {
            hits = results.results.hits;
        } else if (results.hits) {
            hits = results.hits;
        } else if (Array.isArray(results)) {
            hits = results;
        }

        console.log(`Found ${hits.length} transcript hits`);

        // 4. Create summary TSV
        let summaryTSV = 'transcript_id\taccession\tname\tdescription\tevalue\tscore\tbias\tnum_domains\n';

        for (const hit of hits) {
            const name = hit.name || hit.acc || 'N/A';
            const acc = hit.acc || hit.accession || 'N/A';
            const desc = hit.desc || hit.description || 'N/A';
            const evalue = hit.evalue || 'N/A';
            const score = hit.score || 'N/A';
            const bias = hit.bias || 'N/A';
            const numDomains = hit.domains ? hit.domains.length : 0;

            summaryTSV += `${name}\t${acc}\t${name}\t${desc}\t${evalue}\t${score}\t${bias}\t${numDomains}\n`;
        }

        downloadText(summaryTSV, 'hmmer_summary.tsv');
        console.log('Downloaded: hmmer_summary.tsv');

        await delay(500);

        // 5. Create detailed domains table
        let domainsTSV = 'sequence_name\tsequence_acc\tdomain_name\tdomain_acc\tdomain_desc\t';
        domainsTSV += 'seq_start\tseq_end\thmm_start\thmm_end\tievalue\tcevalue\tscore\n';

        for (const hit of hits) {
            const seqName = hit.name || 'N/A';
            const seqAcc = hit.acc || hit.accession || 'N/A';

            const domains = hit.domains || [];
            for (const dom of domains) {
                const domName = dom.alihmmacc || dom.alihmmname || dom.name || 'N/A';
                const domAcc = dom.alihmmname || dom.alihmmacc || dom.acc || 'N/A';
                const domDesc = dom.alihmmdesc || dom.desc || 'N/A';
                const seqFrom = dom.alisqfrom || dom.seq_from || dom.ali_from || 'N/A';
                const seqTo = dom.alisqto || dom.seq_to || dom.ali_to || 'N/A';
                const hmmFrom = dom.alihmmfrom || dom.hmm_from || 'N/A';
                const hmmTo = dom.alihmmto || dom.hmm_to || 'N/A';
                const ievalue = dom.ievalue || dom.i_evalue || 'N/A';
                const cevalue = dom.cevalue || dom.c_evalue || 'N/A';
                const domScore = dom.bitscore || dom.score || 'N/A';

                domainsTSV += `${seqName}\t${seqAcc}\t${domName}\t${domAcc}\t${domDesc}\t`;
                domainsTSV += `${seqFrom}\t${seqTo}\t${hmmFrom}\t${hmmTo}\t${ievalue}\t${cevalue}\t${domScore}\n`;
            }
        }

        downloadText(domainsTSV, 'hmmer_all_domains.tsv');
        console.log('Downloaded: hmmer_all_domains.tsv');

        console.log('');
        console.log('='.repeat(60));
        console.log('DOWNLOAD COMPLETE!');
        console.log('='.repeat(60));
        console.log(`Total transcripts: ${hits.length}`);
        console.log('Files downloaded:');
        console.log('  - hmmer_all_results.json (complete raw data)');
        console.log('  - hmmer_domain_architectures.json (if available)');
        console.log('  - hmmer_summary.tsv (transcript summary)');
        console.log('  - hmmer_all_domains.tsv (all domain predictions)');

    } catch (error) {
        console.error('Error:', error);
        console.log('');
        console.log('If the API endpoint changed, try inspecting network traffic');
        console.log('in the Network tab while clicking around the results page.');
    }
})();
